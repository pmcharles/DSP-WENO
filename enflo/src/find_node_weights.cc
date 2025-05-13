#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include "fv.h"


using namespace std;


//------------------------------------------------------------------------------
// Find Node Weights
//------------------------------------------------------------------------------
void FV::find_node_weights()
{
    
    // Here (i,j,k) refer to node (i-h/2,j-h/2,k-h/2)
    // Setting initial node weights to 1
    for(int i = 0; i< Nx_cell+1; ++i)
        for(int j = 0; j< Ny_cell+1; ++j)
            for(int k = 0; k< Nz_cell+1; ++k)
			    node_w[i][j][k] = 1.0;          

    // If a wall is periodic, divide weight by 2. This can happen seveal times for a single node
    
	// X-Left Wall 
	if(param.bc[0] == periodic)
	   for(int j = 0; j< Ny_cell+1; ++j)
          for(int k = 0; k< Nz_cell+1; ++k)
			    node_w[0][j][k] /= 2.0;	
    
    // X-Right Wall 
	if(param.bc[1] == periodic)
	   for(int j = 0; j< Ny_cell+1; ++j)
          for(int k = 0; k< Nz_cell+1; ++k)
			    node_w[Nx_cell][j][k] /= 2.0;
                
    // Y-Left Wall 
	if(param.bc[2] == periodic)
	   for(int i = 0; i< Nx_cell+1; ++i)
          for(int k = 0; k< Nz_cell+1; ++k)
			    node_w[i][0][k] /= 2.0;	
    
    // Y-Right Wall 
	if(param.bc[3] == periodic)
	   for(int i = 0; i< Nx_cell+1; ++i)
          for(int k = 0; k< Nz_cell+1; ++k)
			    node_w[i][Ny_cell][k] /= 2.0;            	            

	// Z-Left Wall 
	if(param.bc[4] == periodic)
	   for(int i = 0; i< Nx_cell+1; ++i)
          for(int j = 0; j< Ny_cell+1; ++j)
			    node_w[i][j][0] /= 2.0;	
    
    // Z-Right Wall 
	if(param.bc[5] == periodic)
	   for(int i = 0; i< Nx_cell+1; ++i)
          for(int j = 0; j< Ny_cell+1; ++j)
			    node_w[i][j][Nz_cell] /= 2.0; 
}

