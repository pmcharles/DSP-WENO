# Created by: Deep Ray, University of Maryland
# Date: March 25, 2023
# Modified by: Phil Charles (4/15/23)

import os,shutil
import numpy as np
import matplotlib.pyplot as plt

def make_save_dir(PARAMS):
    ''' This function creates the results save directory'''
   
    savedir = get_save_dir(PARAMS)

    if not os.path.exists(savedir):
        os.makedirs(savedir)

    #ckpt_dir = f"{savedir}/ckpts"
    #if not os.path.exists(ckpt_dir):
    #    os.makedirs(ckpt_dir)

    print('\n --- Saving parameters to file \n')
    param_file = savedir + '/parameters.txt'
    with open(param_file,"w") as fid:
        for pname in vars(PARAMS):
            fid.write(f"{pname} = {vars(PARAMS)[pname]}\n")    

    return savedir

def get_save_dir(PARAMS):
    ''' This function create the directory name directory'''

    base_dir  = PARAMS.base_dir
    
    # To ensure no overwriting can occur
    count = 0
    d_name    = f"DSP-WENO_Epoch{PARAMS.max_epoch}_bs{PARAMS.batch_size}_#{count}"
    savedir = f"{base_dir}/{d_name}"
    while os.path.exists(savedir):
        d_name    = f"DSP-WENO_Epoch{PARAMS.max_epoch}_bs{PARAMS.batch_size}_#{count}"
        savedir = f"{base_dir}/{d_name}"
        count += 1
    
    #savedir += PARAMS.dir_suffix

    return savedir     


"""
helper function for plotting
"""
def plot_subplots(f_fullfield, k_fullfield, out_fullfield, u_pred, du_pred, error_vals,savedir,t, sample_type='best'):
    n_samples = f_fullfield.shape[0]
    n_rows = n_samples
    n_cols = 7
    img_size = 32
    NN = img_size**2
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*8,n_rows*8))
    axs = axs.flatten()
    axs_ind = 0
    for n in range(n_samples):
    
        im = axs[axs_ind].imshow(f_fullfield[n,:,:], cmap="inferno")
        fig.colorbar(im, ax=axs[axs_ind],fraction=0.046, pad=0.04)
        axs_ind+=1

        im = axs[axs_ind].imshow(k_fullfield[n,:,:], cmap="inferno")
        fig.colorbar(im, ax=axs[axs_ind],fraction=0.046, pad=0.04)
        axs_ind+=1

        im = axs[axs_ind].imshow(out_fullfield[n,0,:,:], cmap="inferno")
        fig.colorbar(im, ax=axs[axs_ind],fraction=0.046, pad=0.04)
        axs_ind+=1

        im = axs[axs_ind].imshow(u_pred[n,:,:], cmap="inferno")
        fig.colorbar(im, ax=axs[axs_ind],fraction=0.046, pad=0.04)
        axs_ind+=1

        im = axs[axs_ind].imshow(np.abs(u_pred[n,:,:] - out_fullfield[n,0,:,:]), cmap="inferno")
        fig.colorbar(im, ax=axs[axs_ind],fraction=0.046, pad=0.04)
        axs[axs_ind].set_title(f"{error_vals[0][n]:.2e}, {error_vals[1][n]:.2e}")
        axs_ind+=1

        im = axs[axs_ind].imshow(du_pred[n,0,:,:] - out_fullfield[n,1,:,:], cmap="inferno")
        fig.colorbar(im, ax=axs[axs_ind],fraction=0.046, pad=0.04)
        axs_ind+=1

        im = axs[axs_ind].imshow(du_pred[n,1,:,:] - out_fullfield[n,2,:,:], cmap="inferno")
        fig.colorbar(im, ax=axs[axs_ind],fraction=0.046, pad=0.04)
        axs_ind+=1

    fig.savefig(f"{savedir}/{sample_type}_{t+1}.pdf",bbox_inches="tight")
    plt.close()