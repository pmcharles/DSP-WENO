# This is SP-WENO neural network with the sign property strongly enforced

import copy
import time
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as opt
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader
import random
from matplotlib import pyplot as plt
import os, sys
import seaborn as sns
from config import cla
from utils import *
import time

"""
The multi-layer perceptron
"""
class MLP(nn.Module):
    
    def __init__(self, first_layer_input_size, h_widths, act_func, act_func_par):
        super().__init__()
        
        self.first_layer_input_size = first_layer_input_size
        
        # For the output layer
        self.relu = nn.ReLU()
        self.sig  = nn.Sigmoid()
        
        # Activation function
        if act_func == 'ReLU':
            self.act = nn.ReLU()
        elif act_func == 'LeakyReLU':
            self.act = nn.LeakyReLU(act_func_par)
        elif act_func == 'Sigmoid':
            self.act = nn.Sigmoid()
        elif act_func == 'Tanhshrink':
            self.act = nn.Tanhshrink()
        elif act_func == 'ReLU_squared':
            self.act = lambda x: torch.pow(self.relu(x),2)
        else:
            raise Exception("Not a supported activation function.")
        
        
        # Layer sizes (output is 5-dimensional for vertex algorithm)
        self.n_neurons = [self.first_layer_input_size] + h_widths + [5]
        self.n_layers  = len(h_widths) + 1
        
        # Define network architecture
        self.layers = nn.ModuleList()
        for i in range(0,self.n_layers):
            self.layers.append(nn.Linear(self.n_neurons[i], self.n_neurons[i+1])) 

        # Initialize parameters (maybe try standard)
        for i in range(0,self.n_layers):
            nn.init.xavier_uniform_(self.layers[i].weight)
        
    
    def forward(self,state):
        alpha   = 1 - state[:,1] # theta plus is the second data point
        beta    = 1 - state[:,0] # theta minus is the first data point
        
        x_vertices = state[:,5:10]
        y_vertices = state[:,10:15]
        state      = state[:,0:self.first_layer_input_size]
        
        for i in range(self.n_layers):
            state = self.layers[i](state)
            if i < self.n_layers-1:
                state = self.act(state)
        
        # Apply softmax to ensure convex weights
        soft  = nn.Softmax(dim=1)
        state = soft(state)
        
        # Then dot the convex weights with the vertices
        C1 = torch.sum(state*x_vertices,1)
        C2 = torch.sum(state*y_vertices,1)
        
        omega1    = 1/4 - 2*C1
        Tilomega0 = 1/4 - 2*C2
        
        omega1    = torch.transpose(omega1.unsqueeze(0),0,1)
        Tilomega0 = torch.transpose(Tilomega0.unsqueeze(0),0,1)
        
        weights = torch.cat((1-omega1, omega1, Tilomega0, 1-Tilomega0),1)
        return weights, state

# Calculates the reconstructed interface values with the network weights
def inter_evaluation(data_size, sten, cell_centers, weights):
    eno_stencil_values = torch.zeros(data_size,4)
    
    # First two are for the left reconstructed and the other two are for the right
    eno_stencil_values[:,0] = torch.sum(sten[1][:]*cell_centers[:,1:3],dim=1)
    eno_stencil_values[:,1] = torch.sum(sten[2][:]*cell_centers[:,0:2],dim=1)
    eno_stencil_values[:,2] = torch.sum(sten[0][:]*cell_centers[:,2:4],dim=1)
    eno_stencil_values[:,3] = torch.sum(sten[1][:]*cell_centers[:,1:3],dim=1)
    
    network_value_left   = torch.sum(weights[:,:2]*eno_stencil_values[:,:2],dim=1)
    network_value_right  = torch.sum(weights[:,2:]*eno_stencil_values[:,2:],dim=1)
    
    return network_value_left, network_value_right

# Count how many instances where the sign property is satisfied
# original_jump and recon_jump are tensors of the same size
def sign_property(original_jump, recon_jump, data):
    
    cond1 = (original_jump*recon_jump > 0).long()
    cond2 = (torch.abs(original_jump) < 1e-14).long()
    cond3 = (torch.abs(recon_jump) < 1e-14).long()
    
    return sum((torch.logical_or(torch.logical_or(cond1,cond2),cond3)).long())

# This gives a final report on sign property violations from a trained network
# Generates a separate text file for the report
def generate_sign_property_report(original_jump, recon_jump, data, savedir, type_data):
    
    file_name  = f"{savedir}/sign_property_report.txt"
    
    violations = []
    
    # Go through each jump one by one
    for i in range(len(data)):
        # Testing if violation occurs
        #########################################
        acond1 = (original_jump[i]*recon_jump[i] < 0)
        acond2 = (torch.abs(original_jump[i]) > 1e-14)
        acond3 = (torch.abs(recon_jump[i]) > 1e-14)
        aconds = torch.logical_and(torch.logical_and(acond1,acond2),acond3).item()
        if aconds:
            violations.append(i)
        ##########################################
    
    if type_data == 0:
        opening = 'w'
    else:
        opening = 'a'
    
    with open(file_name, opening) as f:
        if type_data == 0:
            f.write('Training Set\n')
            f.write('------------\n\n')
        elif type_data == 1:
            f.write('Validation Set\n')
            f.write('--------------\n\n')
        elif type_data == 2:
            f.write('Test Set\n')
            f.write('--------\n\n')
        
        f.write(f"Number of violation instances: {len(violations)}\n")
        f.write(f"Overall sign property rate: {100*(len(data)-len(violations))/(len(data)):.2f}%\n\n")
        
        if len(violations) > 0:
            f.write('Instances')
            f.write('---------\n\n')
            for point in violations:
                f.write(f"Original jump: {original_jump[point]}\n")
                f.write(f"Reconstructed jump: {recon_jump[point]}\nData:\n")
                np.savetxt(f,data[point,:])
                f.write('\n\n')
            
    return

# Make sure to keep the original values for reconstruction purposes
# This is not neccessary since the scaling is done in the pre-processing of the data in Matlab
def scale_data(data, method):
    scaled_data = data.clone()
    maxes, indices = torch.max(abs(data[:,1:5]),dim=1)
    if method == 1:
        scaling_constants = torch.maximum(torch.ones_like(maxes),maxes)
        scaled_data[:,1:5] /= scaling_constants.view(scaling_constants.size(dim=0),1)
        scaled_data[:,7:10] /= scaling_constants.view(scaling_constants.size(dim=0),1)
        print(torch.max(scaled_data[:,7:10]))
    elif method == 2:
        scaled_data[:,1:5] /= maxes.view(maxes.size(dim=0),1)
        scaled_data[:,7:10] /= maxes.view(maxes.size(dim=0),1)
    elif method == 3:
        maxes, indices = torch.max(data[:,1:5],dim=1)
        scaled_data[:,1:5] /= maxes.view(maxes.size(dim=0),1)
        scaled_data[:,7:10] /= maxes.view(maxes.size(dim=0),1)
    elif method == 4:
        means = torch.mean(scaled_data[:,1:5],1)
        print(means.view(means.size(dim=0),1))
        print(scaled_data[:,1:5])
        scaled_data[:,1:5] -= means.view(means.size(dim=0),1)
        print(scaled_data[:,1:5])
        maxes, indices = torch.max(data[:,1:5],dim=1)
        scaled_data[:,1:5] /= maxes.view(maxes.size(dim=0),1)
        scaled_data[:,7]   /= maxes
    
    return scaled_data

def save_net_weights(mlp, file_name):
    with open(file_name, 'w') as f:
        for i in range(2):
            if i==0:
                f.write('Weights\n')
                for i in range(mlp.n_layers):
                    f.write(f'Layer {i+1}\n')
                    f.write('BEGIN\n')
                    np.savetxt(f,mlp.layers[i].weight.detach().numpy(),delimiter=' ',newline='\n')
                    f.write('END\n')
            else:
                f.write('Biases\n')
                for i in range(mlp.n_layers):
                    f.write(f'Layer {i+1}\n')
                    f.write('BEGIN\n')
                    np.savetxt(f,mlp.layers[i].bias.detach().numpy(),delimiter=' ',newline='\n')
                    f.write('END\n')

def plot_SPWENO_weights(data_points, weights, savedir, name):
    avg_weights  = np.mean(weights,0)
    
    fig = plt.subplots(2, 2, figsize=(24, 24))

    plt.subplot(2, 2, 1)
    plt.plot(data_points, weights[:,0], color = 'green')
    plt.plot(data_points, torch.ones(len(data_points))*avg_weights[0], color = 'black')
    plt.ylabel("Weight value")
    plt.xlabel("Test point")
    plt.legend(['Raw','Average'])
    plt.title("Weight w0 from Network")

    plt.subplot(2, 2, 2)
    plt.plot(data_points, weights[:,1], color = 'green')
    plt.plot(data_points, torch.ones(len(data_points))*avg_weights[1], color = 'black')
    plt.ylabel("Weight value")
    plt.xlabel("Test point")
    plt.legend(['Raw','Average'])
    plt.title("Weight w1 from Network")

    plt.subplot(2, 2, 3)
    plt.plot(data_points, weights[:,2], color = 'red')
    plt.plot(data_points, torch.ones(len(data_points))*avg_weights[2], color = 'black')
    plt.ylabel("Weight value")
    plt.xlabel("Test point")
    plt.legend(['Raw','Average'])
    plt.title("Weight w0-tilde from Network")

    plt.subplot(2, 2, 4)
    plt.plot(data_points, weights[:,3], color = 'red')
    plt.plot(data_points, torch.ones(len(data_points))*avg_weights[3], color = 'black')
    plt.ylabel("Weight value")
    plt.xlabel("Test point")
    plt.legend(['Raw','Average'])
    plt.title("Weight w1-tilde from Network")

    plt.savefig(f"{savedir}/{name}.pdf")
    
torch.set_default_dtype(torch.float64)

#--------------------------------------------------------------
# Loading and saving run parameters
PARAMS = cla()

BATCH_SIZE    = PARAMS.batch_size
EPOCHS        = PARAMS.max_epoch
h_widths      = PARAMS.h_widths
act_func      = PARAMS.act_func
act_func_par  = PARAMS.act_func_par
LR            = PARAMS.lr
loss_func     = PARAMS.loss_func
adam_beta1    = PARAMS.adam_beta1
adam_beta2    = PARAMS.adam_beta2
weight_dec    = PARAMS.weight_dec
val_freq      = PARAMS.val_freq

p_train       = PARAMS.n_train
p_val         = PARAMS.n_val
p_test        = PARAMS.n_test
scaling       = PARAMS.scaling
scale_method  = PARAMS.scale_met

file_name     = PARAMS.data_file
tanh_theta    = PARAMS.tanh_theta

data_pt_size  = 9 + 2 + 10

if tanh_theta:
    data_pt_size += 2
    start_feat = 7
else:
    start_feat = 5
input_size    = PARAMS.input_size 

first_layer_input_size = input_size - 10

print('\n --- Creating operator folder \n')
savedir = make_save_dir(PARAMS)
#--------------------------------------------------------------

torch.manual_seed(PARAMS.seed_no)

#--------------------------------------------------------------
# Load data 

# Unpack data from file
with open(file_name) as f:
    lines = f.readlines()
    text = "".join(lines)

read_data    = []

# Now also labeling the data points for easy access in the raw data file
count = 1
for line in lines:
    tokens  = [t for t in line.split(',')]
    
    temp    = np.zeros(data_pt_size+1+4) # +1 for tag of what data point this is
    
    temp[0] = count
    
    for i in range(1,data_pt_size+1):
        temp[i] = tokens[i-1]
        if i < 5:
            temp[data_pt_size+i] = tokens[i-1]
    
    count  += 1
    
    read_data.append(torch.tensor(temp,dtype=torch.float64))

N = len(read_data)

# Add 4 data points for the unscaled cell center values
data = torch.zeros(N,data_pt_size+1+4,dtype=torch.float64)

for i in range(N):
    data[i,:] = read_data[i]

print(N)

# Scale data if applicable
if scaling:
    scaled_data = scale_data(data,scale_method)
else:
    scaled_data = data.clone()

training_data = scaled_data[:int(p_train*N)]
val_data      = scaled_data[int(p_train*N):int((p_train+p_val)*N),0:data_pt_size-1]
val_actual    = scaled_data[int(p_train*N):int((p_train+p_val)*N),data_pt_size-1:data_pt_size+1]
val_origvals  = scaled_data[int(p_train*N):int((p_train+p_val)*N),data_pt_size+1:data_pt_size+5]
test_data     = scaled_data[int((1-p_test)*N):]

train_dataloader = DataLoader(training_data, batch_size=BATCH_SIZE, shuffle=True)
#----------------------------------------------------------------

# Define the MLP and optimizer
mlp = MLP(first_layer_input_size, h_widths, act_func, act_func_par)
optimizer  = opt.Adam(mlp.parameters(),lr=LR,betas=(adam_beta1,adam_beta2),weight_decay=weight_dec)

if loss_func == 'MSE':
    lossFunc = nn.MSELoss()
elif loss_func == 'L1':
    lossFunc = nn.L1Loss(reduction='mean')
elif loss_func == 'Linf':
    lossFunc = lambda x,y: torch.max(torch.abs(x-y))
elif loss_func == 'MRAE':
    eps = 1e-16
    lossFunc = lambda x,y: torch.mean(torch.abs(x-y)/(torch.abs(y+eps)))

# Train network
mlp.train()

# Store loss for plotting
Loss_epochs = np.zeros(EPOCHS)
Val_loss    = np.zeros(int(EPOCHS/val_freq))

# Store sign property percentage for plotting
Sign_epochs = np.zeros(EPOCHS)
Val_sign    = np.zeros(int(EPOCHS/val_freq))

# Store reconstruction error for plotting
Err_epochs = np.zeros(EPOCHS)
Val_err    = np.zeros(int(EPOCHS/val_freq))

# Stencil
sten = torch.Tensor([[3/2,-1/2],[1/2,1/2],[-1/2,3/2]])

start_time = time.time()
for epoch in range(EPOCHS):
    trainLOSS = 0
    trainERR  = 0
    samples   = 0
    batches   = 0
    
    n_sign_props = 0 # Stores the number of times the sign property is satisfied
    
    for idx, batch in enumerate(train_dataloader):
        data_set = batch[:,0:data_pt_size-1]
        actual   = batch[:,data_pt_size-1:data_pt_size+1]
        orig_val = batch[:,data_pt_size+1:data_pt_size+5] # Unscaled cell center values
        
        weights, v_weights = mlp(data_set[:,start_feat:start_feat+input_size])
        
        network_value_left, network_value_right = inter_evaluation(BATCH_SIZE, sten, orig_val, weights)
        
        n_sign_props += sign_property(orig_val[:,2]-orig_val[:,1],network_value_right-network_value_left,data_set)
        
        # Compare both the left and right reconstructed interface values with the actual value
        err = (lossFunc(network_value_left,actual[:,0]) + lossFunc(network_value_right,actual[:,1]))/2
        
        loss = err
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if loss_func == 'Linf':
            trainLOSS = max(loss.item(),trainLOSS)
            trainERR  = max(err.item(),trainERR)
        else:
            trainLOSS += loss.item()
            trainERR  += err.item()
        batches   += 1
        samples   += BATCH_SIZE
    
    # Use the validation set
    if (epoch+1) % val_freq == 0:
        mlp.eval()
        weights, v_weights = mlp(val_data[:,start_feat:start_feat+input_size])
        
        network_value_left, network_value_right = inter_evaluation(len(val_data), sten, val_origvals, weights)
        
        n_val_sign_props = sign_property(val_origvals[:,2]-val_origvals[:,1],network_value_right-network_value_left, val_data)
        
        # Compare both the left and right reconstructed interface values with the actual value
        valERR = (lossFunc(network_value_left,val_actual[:,0]) + lossFunc(network_value_right,val_actual[:,1]))/2 
        
        valLOSS = valERR
        
        print("Epoch: {}\t Validation Loss: {}\t Validation Error: {}\t Percent Sign Property: {}".format(
                                                      epoch+1, valLOSS, valERR, (n_val_sign_props/len(val_data))))
        
        Val_loss[int((epoch+1)/val_freq)-1] = valLOSS
        Val_err[int((epoch+1)/val_freq)-1]  = valERR
        Val_sign[int((epoch+1)/val_freq)-1] = n_val_sign_props/len(val_data)
        
        mlp.train()
    
    if loss_func == 'Linf':
        Loss_epochs[epoch] = trainLOSS
        Err_epochs[epoch]  = trainERR
    else:
        Loss_epochs[epoch] = trainLOSS/batches
        Err_epochs[epoch]  = trainERR/batches
    Sign_epochs[epoch] = n_sign_props/samples
    
end_time = time.time()

# Evaluate network
mlp.eval()

# Save network in PyTorch
name = "trained_network"
PATH = f"{savedir}/{name}"
torch.save(mlp.state_dict(), PATH)

# Save final results for training, validation, and test data after training
file_name = f"{savedir}/after_training_results.txt"

with open(file_name, 'w') as f:
    for i in range(3):
        
        if i == 0:
            data_set = training_data[:,0:data_pt_size-1]
            actual   = training_data[:,data_pt_size-1:data_pt_size+1]
            orig_val = training_data[:,data_pt_size+1:data_pt_size+5] # Unscaled cell center values
            
            f.write('Training Set\n')
            f.write('--------------\n')
        elif i == 1:
            data_set = val_data
            actual   = val_actual
            orig_val = val_origvals
            
            f.write('Validation Set\n')
            f.write('--------------\n')
        elif i == 2:
            data_set = test_data[:,0:data_pt_size-1]
            actual   = test_data[:,data_pt_size-1:data_pt_size+1]
            orig_val = test_data[:,data_pt_size+1:data_pt_size+5]
            
            f.write('Test Set\n')
            f.write('--------------\n')
        
        weights, conv_weights = mlp(data_set[:,start_feat:start_feat+input_size])
        
        network_value_left, network_value_right = inter_evaluation(len(data_set), sten, orig_val, weights)
        left_error  = lossFunc(network_value_left,actual[:,0])
        right_error = lossFunc(network_value_right,actual[:,1])
        perc_sign_props = sign_property(orig_val[:,2]-orig_val[:,1],network_value_right-network_value_left, data_set)/len(data_set) * 100
        
        generate_sign_property_report(orig_val[:,2]-orig_val[:,1],network_value_right-network_value_left, data_set, savedir, i)
        
        f.write(f"Left Reconstruction Error: {left_error}\n")
        f.write(f"Right Reconstruction Error: {right_error}\n")
        f.write(f"Sign Property: {perc_sign_props:.2f}%\n\n")
    
    f.write(f"Training time: {end_time-start_time} s")

# Save network to be accessed in Matlab
name = "dsp_weno"
file_name = f"{savedir}/{name}.txt"
save_net_weights(mlp, file_name)

# Plotting results
#################################################
sns.set()
sns.set_style("dark")

data_points  = torch.arange(len(test_data))
test_weights = weights.detach().numpy()
avg_weights  = torch.mean(weights.detach(),0)
#print(avg_weights)
epochs   = np.arange(1,EPOCHS+1,1)
v_epochs = np.arange(val_freq, EPOCHS+1, val_freq)

fig = plt.subplots(2, 2, figsize=(24, 24))

plt.subplot(2, 2, 1)
plt.semilogy(epochs,Err_epochs, color = 'navy')
plt.semilogy(v_epochs, Val_err, color = 'orange')
plt.title('Network Error (Log Scale)')
plt.legend(['Training','Validation'])
plt.xlabel("Epochs")
plt.ylabel("Error")

plt.subplot(2, 2, 2)
plt.semilogy(epochs, Loss_epochs, color = 'navy')
plt.semilogy(v_epochs, Val_loss, color = 'orange')
plt.title('Network Log Loss')
plt.legend(['Training','Validation'])
plt.xlabel("Epochs")
plt.ylabel("Log Loss")

plt.subplot(2, 2, 3)
plt.plot(epochs, Sign_epochs, color = 'navy')
plt.title('Sign Property (Training)')
plt.xlabel("Epochs")
plt.ylabel("Sign Property Percentage")

plt.subplot(2, 2, 4)
plt.plot(v_epochs, Val_sign, color = 'navy')
plt.title('Sign Property (Validation)')
plt.xlabel("Epochs")
plt.ylabel("Sign Property Percentage")

plt.savefig(f"{savedir}/results.pdf")


# Plot the SP-WENO weights generated from the NN
#################################################
plot_SPWENO_weights(data_points, weights.detach().numpy(), savedir, 'SP-WENO_weights')