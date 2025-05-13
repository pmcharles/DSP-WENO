#include <vector>
#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <cmath>
#include <algorithm>
#include "neural_network.h"

using namespace std;

//------------------------------------------------------------------------------
// The input is assumed to be just the cell center values. All scaling and
// computation of derived quantities is done here
//------------------------------------------------------------------------------
void Neural_Network::forward(const double input[], double &recl, double &recr) const
{
    // Perhaps move this into the param file
    int vertex_mode  = 2;

    double cell_cen_vals[4];
    for (int i = 0; i < 4; i++)
        cell_cen_vals[i] = input[i];

    scale(cell_cen_vals);

    const double theta_m = (cell_cen_vals[3]-cell_cen_vals[2])/(cell_cen_vals[2]-cell_cen_vals[1]);
    const double theta_p = (cell_cen_vals[1]-cell_cen_vals[0])/(cell_cen_vals[2]-cell_cen_vals[1]);

    double del1 = abs(cell_cen_vals[1]-cell_cen_vals[0]);
    double del2 = abs(cell_cen_vals[2]-cell_cen_vals[1]);
    double del3 = abs(cell_cen_vals[3]-cell_cen_vals[2]);

    double max_del = max(max(del1,del2),del3);
    
    vector<double> state = { tanh(theta_m), tanh(theta_p), del1, del2, del3 };
    vector< vector<double> > vertices = compute_vertices(theta_m,theta_p,max_del,vertex_mode);

    // Now run the neural network
    for (int i = 0; i < n_layers; i++)
    {
        vector<double> temp;
        for (int j = 0; j < network_weights[i].size(); j++)
        {
            double dot = 0.0;
            for (int k = 0; k < network_weights[i][j].size(); k++)
            {
                dot += network_weights[i][j][k]*state[k];
            }
            temp.push_back(dot);
        }

        for (int j = 0; j < network_biases[i].size(); j++)
        {
            state[j] = temp[j] + network_biases[i][j];
        }

        if (i != n_layers-1)
            relu(state);
    }

    softmax(state);
    // Computing the perturbations
    double C1 = 0.0, C2 = 0.0;
    for (int i = 0; i < state.size(); i++)
    {
        C1 += state[i]*vertices[0][i];
        C2 += state[i]*vertices[1][i];
    }

    double w1      = 1.0/4.0 - 2*C1;
    double Tildew0 = 1.0/4.0 - 2*C2;

    vector<double> weights = { 1.0-w1, w1, Tildew0, 1.0-Tildew0 };

    compute_recon(weights, input, recl, recr);
}

//------------------------------------------------------------------------------
// Load Neural Network weights and biases from txt file
//------------------------------------------------------------------------------
void Neural_Network::load_network(const std::string net_file_name)
{
	ifstream in(net_file_name.c_str());
	if (in.good() == 0) {
        cout << "DSP-WENO network file does not exist! Exiting!" << endl;
        exit(0);
    }
    
	int weight = 1;
	string temp;
	while (getline(in,temp))
	{
		if (temp.compare("BEGIN") == 0)
		{
			vector< vector<double> > temp_mat;
			vector<double> temp_bias;
			
			getline(in,temp);
			while (temp.compare("END") != 0)
			{
				vector<double> temp_row;
				istringstream iss(temp);

				for (;;)
				{
					double val;
					iss >> val;
					if (!iss) break;

					if (weight == 1)
						temp_row.push_back(val);
					else
						temp_bias.push_back(val);
				}
				if (weight == 1)
					temp_mat.push_back(temp_row);
				getline(in,temp);
			}
			if (weight == 1)
				network_weights.push_back(temp_mat);
			else
				network_biases.push_back(temp_bias);
		}
		else if (temp.compare("Biases") == 0)
		{
			weight = 0;
		}
	}
}

//------------------------------------------------------------------------------
// Computes the vertices of the applicable convex shape in the feasible region
//------------------------------------------------------------------------------
std::vector< std::vector<double> > Neural_Network::compute_vertices(const double theta_m, const double theta_p, const double max_delta, const int vertex_mode, const int permutation) const
{
    double threshold = 1e-8;

    double psi_p = (1-theta_m)/(1-theta_p);
    double psi_m = 1/psi_p;

    // To ensure sum_deltas doesn't exceed the feasible region
    double gamma1 =  min(max_delta,1.0/8.0);
    double gamma2 = -min(max_delta,3.0/8.0);
    
    // x1 & y1 correspond to when psi_p < -1 and x2 & y2 otherwise in
    // cases 2 and 3
    double x1 = 1.0/8.0*(1+psi_p) - gamma1*psi_p;
    double y1 = 1.0/8.0*(1+psi_m) - gamma2*psi_m;
    double x2 = 1.0/8.0*(1+psi_p) - gamma2*psi_p;
    double y2 = 1.0/8.0*(1+psi_m) - gamma1*psi_m;

    vector< vector<double> > vertices;

    if ((theta_m - 1.0) > threshold && (theta_p-1.0) > threshold)
    {
        vertices = { {1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0}, {1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0} };
    }
    else if ((theta_m - 1.0) > threshold && (theta_p - 1.0) < -1.0*threshold)
    {
        if (psi_p < -1.0)
        {            
            if (x1 < gamma2)
            {
                // Just a single point (closest point to origin)
                double new_x = ((psi_m*psi_m)+psi_m)/(8*((psi_m*psi_m)+1));
                double new_y =  1.0/8.0*(1+psi_m) - new_x*psi_m;
                vertices = { {new_x, new_x, new_x, new_x, new_x}   , {new_y, new_y, new_y, new_y, new_y} };
            }
            else
            {
                if (vertex_mode == 1 || vertex_mode == 2) // old
                {
                    double cent_x = (2*gamma2 + x1)/3.0;
                    double cent_y = (2*gamma1 + y2)/3.0;
                    vertices = { {gamma2, x1, gamma2, cent_x, cent_x}, {gamma1, gamma1, y2, cent_y, cent_y} };
                }
                else if (vertex_mode == 3)
                {
                    double cent_x = (2*gamma2 + max(-3.0/8.0,2*gamma2-x1))/3.0;
                    double cent_y = (2*gamma1 + y2)/3.0;
                    vertices = { {gamma2, max(-3.0/8.0,2*gamma2-x1), gamma2, cent_x, cent_x}, {gamma1, gamma1, y2, cent_y, cent_y} };
                }
            }
        }
        else
        {
            if (y2 < gamma2)
                // Entire O(h) region is good
                vertices = { {gamma2, gamma1, gamma2, gamma1, 0.0}, {gamma1, gamma1, gamma2, gamma2, 0.0} };
            else
                vertices = { {gamma2, gamma1, gamma2, gamma1, x2} , {gamma1, gamma1, gamma2, y2, gamma2} };
        }
    }
    else if ((theta_m - 1.0) < -1.0*threshold && (theta_p - 1.0) > threshold)
    {
        if (psi_p < -1.0)
        {
            if (x1 < gamma2)
                // Entire O(h) region is good
                vertices = { {gamma2, gamma1, gamma2, gamma1, 0.0}, {gamma1, gamma1, gamma2, gamma2, 0.0} };
            else
                vertices = { {gamma1, gamma1, gamma2, x1, gamma2} , {gamma2, gamma1, gamma2, gamma1, y1} };
        }
        else
        {       
                if (y2 < gamma2)
                {
                    // Just a single point (closest point to origin)
                    double new_x = ((psi_m*psi_m)+psi_m)/(8*((psi_m*psi_m)+1));
                    double new_y =  1.0/8.0*(1+psi_m) - new_x*psi_m;
                    vertices = { {new_x, new_x, new_x, new_x, new_x}   , {new_y, new_y, new_y, new_y, new_y} };
                }
                else
                {
                    if (vertex_mode == 1)
                    {
                        double cent_x = (2*gamma1 + x2)/3.0;
                        double cent_y = (2*gamma2 + y2)/3.0;
                        vertices = { {gamma1, gamma1, x2, cent_x, cent_x}, {gamma2, y2, gamma2, cent_y, cent_y} };
                    }
                    else if (vertex_mode == 2)
                    {
                        double cent_x = (2*gamma1 + x1)/3.0;
                        double cent_y = (2*gamma2 + y2)/3.0;
                        vertices = { {gamma1, x1, gamma1, cent_x, cent_x}, {gamma2, gamma2, y2, cent_y, cent_y} };
                    }
                    else if (vertex_mode == 3)
                    {
                        double cent_x = (2*gamma1 + x1)/3.0;
                        double cent_y = (2*gamma2 + max(-3.0/8.0,2*gamma2-y2))/3.0;
                        vertices = { {gamma1, gamma1, x1, cent_x, cent_x}, {gamma2, max(-3.0/8.0,2*gamma2-y2), gamma2, cent_y, cent_y} };
                    }
                }
        }
    }
    else if (abs(theta_m - 1.0) < threshold && (theta_p - 1.0) > threshold)
    {
        if (permutation == 2)
            vertices = { {1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0,  1.0/8.0}, {-3.0/8.0, 1.0/8.0, 1.0/8.0, -3.0/8.0, -1.0/8.0} };
        else
            vertices = { {1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0,  1.0/8.0}, {-3.0/8.0, 1.0/8.0, -3.0/8.0, 1.0/8.0, -1.0/8.0} };
    }
    else if (abs(theta_p - 1.0) < threshold && (theta_m - 1.0) > threshold)
    {
        if (permutation == 2)
            vertices = { {1.0/8.0, 1.0/8.0, -3.0/8.0, -3.0/8.0, -1.0/8.0}, {1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0} };
        else
            vertices = { {-3.0/8.0, 1.0/8.0, -3.0/8.0, 1.0/8.0, -1.0/8.0}, {1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0} };
    }
    else
    {
        if (permutation == 2)
            vertices = { {1.0/8.0, 1.0/8.0, -3.0/8.0, -3.0/8.0, -1.0/8.0}, {1.0/8.0, -3.0/8.0, -3.0/8.0, 1.0/8.0, -1.0/8.0} };
        else
            vertices = { {-3.0/8.0, 1.0/8.0, 1.0/8.0, -3.0/8.0, -1.0/8.0}, {1.0/8.0, 1.0/8.0, -3.0/8.0, -3.0/8.0, -1.0/8.0} };
    }

    return vertices;
}

//------------------------------------------------------------------------------
// ReLU Activation
//------------------------------------------------------------------------------
void Neural_Network::relu(std::vector<double>& state) const
{
    for (int i = 0; i < state.size(); i++)
    {
        state[i] = max(0.0,state[i]);
    }
}

//------------------------------------------------------------------------------
// Computes Softmax of input
//------------------------------------------------------------------------------
void Neural_Network::softmax(std::vector<double>& state) const
{
    double max_val = 0.0;
    for (int i = 0; i < state.size(); i++)
    {
        max_val = max(max_val,state[i]);
    }

    double sum = 0;
    for (int i = 0; i < state.size(); i++)
    {
        state[i] = exp(state[i]-max_val);
        sum += state[i];
    }

    for (int i = 0; i < state.size(); i++)
    {
        state[i] /= sum;
    }
}

//------------------------------------------------------------------------------
// Scales the cell center values so that they lie in [-1,1]
//------------------------------------------------------------------------------
void Neural_Network::scale(double state[]) const
{
    double max_val  = 1.0;
    for (int i = 0; i < 4; i++)
    {
        max_val = max(max_val,abs(state[i]));
    }


    for (int i = 0; i < 4; i++)
    {
        state[i] /= max_val;
    }
}


//------------------------------------------------------------------------------
// Computes the left and right reconstructions given the output weights from
// the network
//------------------------------------------------------------------------------
void Neural_Network::compute_recon(const std::vector<double> weights, const double state[],
                                                  double &recl, double &recr) const
{
    recl = weights[0]*0.5*(state[1] + state[2]) + weights[1]*(-0.5*state[0] + 1.5*state[1]);
    recr = weights[3]*0.5*(state[1] + state[2]) + weights[2]*(-0.5*state[3] + 1.5*state[2]);
}