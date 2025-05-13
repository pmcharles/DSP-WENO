#ifndef __NEURAL_NETWORK_H__
#define __NEURAL_NETWORK_H__

//------------------------------------------------------------------------------
// Neural Network class
//------------------------------------------------------------------------------

#include <vector>
#include <fstream>
#include <string>
#include "ext_constants.h"

class Neural_Network
{
    public:
        Neural_Network() {};
        Neural_Network(const std::string net_file_name, const ReconstructionScheme reconstruct_scheme)
        {
            if (reconstruct_scheme == sp_weno_dl)
            {
                load_network(net_file_name);
                n_layers = network_weights.size();
            }
        };
        void forward(const double state[], double &recl, double &recr) const;

    private:
        int n_layers;
        std::vector< std::vector< std::vector<double> > > network_weights;
	    std::vector< std::vector<double> > 		          network_biases;
        void load_network(const std::string net_file_name);
        std::vector< std::vector<double> > compute_vertices(const double theta_m, const double theta_p, const double max_delta, const int vertex_mode, const int permutation = 2) const;
        void relu(std::vector<double>& state) const;
        void softmax(std::vector<double>& state) const;
        void scale(double state[]) const;
        void compute_recon(const std::vector<double> weights, const double state[], double &recl, double &recr) const;

};

#endif