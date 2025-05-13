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
#include "neural_network.h"

using namespace std;

//------------------------------------------------------------------------------
// Initializing problem
//---------------------------------------------------------------------------
void FV::initialize ()
{
   // Construct the Neural Network
   net = Neural_Network(param.net_file,param.reconstruct_scheme);

   // Make grid
   cout<<"  ---Initializing mesh ...\n";
   
   ng_cell = max(param.flux_ng_cell,param.rec_ng_cell);
   rec_ng_cell = param.rec_ng_cell;
   
   dx = (param.xmax - param.xmin) / Nx_cell;
   dy = (param.ymax - param.ymin) / Ny_cell;
   dz = (param.zmax - param.zmin) / Nz_cell;
   cell_cc.resize (Nx_cell);
   for(int i=0; i<Nx_cell; ++i)
   {
      cell_cc[i].resize(Ny_cell);
      for(int j=0; j<Ny_cell; ++j)
      {
         cell_cc[i][j].resize(Nz_cell);
         for(int k=0; k<Nz_cell; ++k)
		 {
			cell_cc[i][j][k].x = param.xmin + 0.5*dx + i*dx;
			cell_cc[i][j][k].y = param.ymin + 0.5*dy + j*dy;
			cell_cc[i][j][k].z = param.zmin + 0.5*dz + k*dz;
		 }   
      }
   }
   
   // Initial condition
   cout << "  ---Setting initial condition ...\n";
   primitive.resize (Nx_cell + 2*ng_cell);
   residual.resize (Nx_cell);
   if(param.time_scheme == rk4)
      residual2.resize (Nx_cell);
   conserved.resize (Nx_cell);
   conserved_old.resize (Nx_cell);
   entropy_var.resize (Nx_cell + 2*ng_cell);
   vorticity.resize(Nx_cell+1); //At nodes
   node_w.resize(Nx_cell + 1);
   
   for(int i=0; i<Nx_cell + 2*ng_cell; ++i)
   {
	  primitive[i].resize (Ny_cell + 2*ng_cell, vector<PrimVar>(Nz_cell + 2*ng_cell));
	  entropy_var[i].resize (Ny_cell + 2*ng_cell, vector<EntVar>(Nz_cell + 2*ng_cell));
	  if(i>=ng_cell && i < Nx_cell + ng_cell)
	  {
		 residual[i-ng_cell].resize (Ny_cell, vector<Flux>(Nz_cell));
		 if(param.time_scheme == rk4)
            residual2[i-ng_cell].resize (Ny_cell, vector<Flux>(Nz_cell));
		 conserved[i-ng_cell].resize (Ny_cell, vector<ConVar>(Nz_cell));
		 conserved_old[i-ng_cell].resize (Ny_cell, vector<ConVar>(Nz_cell));
	  }
   }
   
   for(int i=0; i<Nx_cell+1; ++i)
   {
	  vorticity[i].resize (Ny_cell + 1, vector<Vector>(Nz_cell+1));
      node_w[i].resize (Ny_cell + 1, vector<double>(Nz_cell+1));
   }  
   find_node_weights();
	  
   for(int i=0; i<Nx_cell; ++i)
   {
      for(int j=0; j<Ny_cell; ++j)
      {
         for(int k=0; k<Nz_cell; ++k)
		 {
			param.initial_condition.value (cell_cc[i][j][k],primitive[i+ng_cell][j+ng_cell][k+ng_cell]);
			assert (primitive[i+ng_cell][j+ng_cell][k+ng_cell].density > 0.0);
			assert (primitive[i+ng_cell][j+ng_cell][k+ng_cell].pressure > 0.0);
		 }   
      }
   }
   
   Cp  = (param.GAMMA * param.gas_const / (param.GAMMA - 1.0));
   
   apply_bc();
}

//------------------------------------------------------------------------------
// Clear data structured for next refinement 
//------------------------------------------------------------------------------
void FV::clear_data_structs ()
{
   // Make grid
   cout<<"  ---Clearing structures ...\n";
   cell_cc.clear();
   node_w.clear();
   vorticity.clear();
   primitive.clear();
   residual.clear();
   conserved.clear();
   conserved_old.clear();
   entropy_var.clear();
   
}


//------------------------------------------------------------------------------
// Apply Boundary condition
//------------------------------------------------------------------------------
void FV::apply_bc()
{
   
    
    // X-Left Wall
	if(param.bc[0] == periodic)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			    primitive[i][ng_cell + j][ng_cell + k] = primitive[Nx_cell + i][ng_cell + j][ng_cell + k];
	}
	else if(param.bc[0] == wall)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				primitive[i][ng_cell + j][ng_cell + k].density    = primitive[2*ng_cell - 1 - i][ng_cell + j][ng_cell + k].density;	
				primitive[i][ng_cell + j][ng_cell + k].velocity.x = -primitive[2*ng_cell - 1 - i][ng_cell + j][ng_cell + k].velocity.x;
				primitive[i][ng_cell + j][ng_cell + k].velocity.y = primitive[2*ng_cell - 1 - i][ng_cell + j][ng_cell + k].velocity.y;
				primitive[i][ng_cell + j][ng_cell + k].velocity.z = primitive[2*ng_cell - 1 - i][ng_cell + j][ng_cell + k].velocity.z;
				primitive[i][ng_cell + j][ng_cell + k].pressure   = primitive[2*ng_cell - 1 - i][ng_cell + j][ng_cell + k].pressure;
			 }
	}
	else if(param.bc[0] == neumann)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				primitive[i][ng_cell + j][ng_cell + k] = primitive[ng_cell][ng_cell + j][ng_cell + k];
			 }
	}
	
	// X-Right Wall
	if(param.bc[1] == periodic)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			    primitive[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k] = primitive[ng_cell + i][ng_cell + j][ng_cell + k];	
	}
	else if(param.bc[1] == wall)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				primitive[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k].density    = primitive[Nx_cell + ng_cell - 1 - i][ng_cell + j][ng_cell + k].density;	
				primitive[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k].velocity.x = - primitive[Nx_cell + ng_cell - 1 - i][ng_cell + j][ng_cell + k].velocity.x;
				primitive[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k].velocity.y = primitive[Nx_cell + ng_cell - 1 - i][ng_cell + j][ng_cell + k].velocity.y;
				primitive[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k].velocity.z = primitive[Nx_cell + ng_cell - 1 - i][ng_cell + j][ng_cell + k].velocity.z;
				primitive[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k].pressure   = primitive[Nx_cell + ng_cell - 1 - i][ng_cell + j][ng_cell + k].pressure;
			 }
	}
	else if(param.bc[1] == neumann)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				primitive[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k] = primitive[Nx_cell + ng_cell - 1][ng_cell + j][ng_cell + k];
			 }
	}
    
    // Y-Left Wall
	if(param.bc[2] == periodic)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			    primitive[ng_cell + i][j][ng_cell + k] = primitive[ng_cell + i][Ny_cell + j][ng_cell + k];
	}
	else if(param.bc[2] == wall)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				primitive[ng_cell + i][j][ng_cell + k].density    = primitive[ng_cell + i][2*ng_cell - 1 - j][ng_cell + k].density;	
				primitive[ng_cell + i][j][ng_cell + k].velocity.x = primitive[ng_cell + i][2*ng_cell - 1 - j][ng_cell + k].velocity.x;
				primitive[ng_cell + i][j][ng_cell + k].velocity.y = - primitive[ng_cell + i][2*ng_cell - 1 - j][ng_cell + k].velocity.y;
				primitive[ng_cell + i][j][ng_cell + k].velocity.z = primitive[ng_cell + i][2*ng_cell - 1 - j][ng_cell + k].velocity.z;
				primitive[ng_cell + i][j][ng_cell + k].pressure   = primitive[ng_cell + i][2*ng_cell - 1 - j][ng_cell + k].pressure;
			 }
	}
	else if(param.bc[2] == neumann)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				primitive[ng_cell + i][j][ng_cell + k] = primitive[ng_cell + i][ng_cell][ng_cell + k];
			 }		  
	} 
    
    // Y-Right Wall
	if(param.bc[3] == periodic)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
			 for(int k = 0; k< Nz_cell; ++k)
			    primitive[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k] = primitive[ng_cell + i][ng_cell + j][ng_cell + k];
	}
	else if(param.bc[3] == wall)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				primitive[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k].density    = primitive[ng_cell + i][Ny_cell + ng_cell - 1 - j][ng_cell + k].density;	
				primitive[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k].velocity.x = primitive[ng_cell + i][Ny_cell + ng_cell - 1 - j][ng_cell + k].velocity.x;
				primitive[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k].velocity.y = - primitive[ng_cell + i][Ny_cell + ng_cell - 1 - j][ng_cell + k].velocity.y;
				primitive[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k].velocity.z = primitive[ng_cell + i][Ny_cell + ng_cell - 1 - j][ng_cell + k].velocity.z;
				primitive[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k].pressure   = primitive[ng_cell + i][Ny_cell + ng_cell - 1 - j][ng_cell + k].pressure;
			 }
	}
	else if(param.bc[3] == neumann)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				primitive[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k] = primitive[ng_cell + i][Ny_cell + ng_cell - 1][ng_cell + k];
			 }
		  
	   
	}
	
	
	// Z-Left Wall
	if(param.bc[4] == periodic)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			    primitive[ng_cell + i][ng_cell + j][k] = primitive[ng_cell + i][ng_cell + j][Nz_cell + k];
	}
	else if(param.bc[4] == wall)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			 {
				primitive[ng_cell + i][ng_cell + j][k].density    = primitive[ng_cell + i][ng_cell + j][2*ng_cell - 1 - k].density;	
				primitive[ng_cell + i][ng_cell + j][k].velocity.x = primitive[ng_cell + i][ng_cell + j][2*ng_cell - 1 - k].velocity.x;
				primitive[ng_cell + i][ng_cell + j][k].velocity.y = primitive[ng_cell + i][ng_cell + j][2*ng_cell - 1 - k].velocity.y;
				primitive[ng_cell + i][ng_cell + j][k].velocity.z = -primitive[ng_cell + i][ng_cell + j][2*ng_cell - 1 - k].velocity.z;
				primitive[ng_cell + i][ng_cell + j][k].pressure   = primitive[ng_cell + i][ng_cell + j][2*ng_cell - 1 - k].pressure;
			 }
	}
	else if(param.bc[4] == neumann)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			 {
				primitive[ng_cell + i][ng_cell + j][k] = primitive[ng_cell + i][ng_cell + j][ng_cell];
			 }		  
	} 
	
	// Z-Right Wall
	if(param.bc[5] == periodic)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
			 for(int k = 0; k< ng_cell; ++k)
			    primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k] = primitive[ng_cell + i][ng_cell + j][ng_cell + k];
	}
	else if(param.bc[5] == wall)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			 {
				primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k].density    = primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell - 1 - k].density;	
				primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k].velocity.x = primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell - 1 - k].velocity.x;
				primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k].velocity.y = primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell - 1 - k].velocity.y;
				primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k].velocity.z = -primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell - 1 - k].velocity.z;
				primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k].pressure   = primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell - 1 - k].pressure;
			 }
	}
	else if(param.bc[5] == neumann)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			 {
				primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k] = primitive[ng_cell + i][ng_cell + j][Nz_cell + ng_cell - 1];
			 }
		  
	   
	}
    
    
    // Wedging edges also needed for viscous flux
    // For the time being, it is set assumin all bc are periodic Need to fix this.
    for(int i=0;i<Nx_cell; ++i)
    {
         primitive[ng_cell + i][ng_cell - 1][ng_cell - 1]             = primitive[ng_cell + i][ng_cell + Ny_cell - 1][ng_cell + Nz_cell - 1];
         primitive[ng_cell + i][ng_cell + Ny_cell][ng_cell- 1]        = primitive[ng_cell + i][ng_cell][ng_cell + Nz_cell - 1];
         primitive[ng_cell + i][ng_cell - 1][ng_cell + Nz_cell]       = primitive[ng_cell + i][ng_cell + Ny_cell - 1][ng_cell];
         primitive[ng_cell + i][ng_cell + Ny_cell][ng_cell + Nz_cell] = primitive[ng_cell + i][ng_cell][ng_cell];
    }
    for(int j=0;j<Ny_cell; ++j)
    {
         primitive[ng_cell - 1][ng_cell + j][ng_cell - 1]             = primitive[ng_cell + Nx_cell - 1][ng_cell + j][ng_cell + Nz_cell - 1];
         primitive[ng_cell + Nx_cell][ng_cell + j][ng_cell - 1]       = primitive[ng_cell][ng_cell + j][ng_cell + Nz_cell - 1];
         primitive[ng_cell - 1][ng_cell + j][ng_cell + Nz_cell]       = primitive[ng_cell + Nx_cell - 1][ng_cell + j][ng_cell];
         primitive[ng_cell + Nx_cell][ng_cell + j][ng_cell + Nz_cell] = primitive[ng_cell][ng_cell + j][ng_cell];
    } 
    for(int k=0;k<Nz_cell; ++k)
    {
         primitive[ng_cell - 1][ng_cell - 1][ng_cell + k]             = primitive[ng_cell + Nx_cell - 1][ng_cell + Ny_cell - 1][ng_cell + k];
         primitive[ng_cell + Nx_cell][ng_cell - 1][ng_cell + k]       = primitive[ng_cell][ng_cell + Ny_cell - 1][ng_cell + k];
         primitive[ng_cell - 1][ng_cell + Ny_cell][ng_cell + k]       = primitive[ng_cell + Nx_cell - 1][ng_cell][ng_cell + k];
         primitive[ng_cell + Nx_cell][ng_cell + Ny_cell][ng_cell + k] = primitive[ng_cell][ng_cell][ng_cell + k];
    }    
      

    // Corner also needed for viscous flux. 
    // -1,-1,-1
	primitive[ng_cell -1][ng_cell -1][ng_cell -1] = primitive[ng_cell + Nx_cell -1][ng_cell + Ny_cell -1][ng_cell + Nz_cell -1];
    
    primitive[ng_cell + Nx_cell][ng_cell -1][ng_cell -1] = primitive[ng_cell][ng_cell + Ny_cell -1][ng_cell + Nz_cell -1];
    
    primitive[ng_cell -1][ng_cell + Ny_cell][ng_cell -1] = primitive[ng_cell + Nx_cell -1][ng_cell][ng_cell + Nz_cell -1];
    
    primitive[ng_cell -1][ng_cell -1][ng_cell + Nz_cell] = primitive[ng_cell + Nx_cell -1][ng_cell + Ny_cell -1][ng_cell];
    
    primitive[ng_cell + Nx_cell][ng_cell + Ny_cell][ng_cell -1] = primitive[ng_cell][ng_cell][ng_cell + Nz_cell -1];
    
    primitive[ng_cell + Nx_cell][ng_cell -1][ng_cell + Nz_cell] = primitive[ng_cell][ng_cell + Ny_cell -1][ng_cell];
    
    primitive[ng_cell -1][ng_cell + Ny_cell][ng_cell + Nz_cell] = primitive[ng_cell + Nx_cell -1][ng_cell][ng_cell];
    
    primitive[ng_cell + Nx_cell][ng_cell + Ny_cell][ng_cell + Nz_cell] = primitive[ng_cell][ng_cell][ng_cell];
    
	
}

//------------------------------------------------------------------------------
// Compute time step (DX and DY factor???)
//------------------------------------------------------------------------------
void FV::compute_dt ()
{
   dt = 1.0e20;
   for(int i=0; i<Nx_cell; ++i)
      for(int j=0; j<Ny_cell; ++j)
         for(int k=0; k<Nz_cell; ++k)
		 {
         double speed = primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity.norm() + 
			   sqrt(param.GAMMA * primitive[i+ng_cell][j+ng_cell][k+ng_cell].pressure / primitive[i+ng_cell][j+ng_cell][k+ng_cell].density);
			if(param.model == ns_kep || param.model == ns_es)
               speed += param.Viscosity(Temperature(primitive[i+ng_cell][j+ng_cell][k+ng_cell]))/ 
                        primitive[i+ng_cell][j+ng_cell][k+ng_cell].density / min(dx,min(dy,dz));
            dt = min (dt, min(dx,min(dy,dz))/speed);   
		 }
   
   dt *= param.cfl;
   
}

//------------------------------------------------------------------------------
// Convert conserved to primitive
//------------------------------------------------------------------------------
void FV::con_to_prim ()
{
   for(int i=0; i<Nx_cell; ++i)
      for(int j=0; j<Ny_cell; ++j)
         for(int k=0; k<Nz_cell; ++k)
		 {
			primitive[i+ng_cell][j+ng_cell][k+ng_cell].density = conserved[i][j][k].density;
			primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity.equ(conserved[i][j][k].momentum,1.0/conserved[i][j][k].density);
			primitive[i+ng_cell][j+ng_cell][k+ng_cell].pressure = (param.GAMMA-1.0)*(conserved[i][j][k].energy 
									   - 0.5*conserved[i][j][k].momentum.square()/conserved[i][j][k].density);
		 }
   apply_bc();
}

//------------------------------------------------------------------------------
// Convert primitive to conserved
//------------------------------------------------------------------------------
void FV::prim_to_con ()
{
   for(int i=0; i<Nx_cell; ++i)
      for(int j=0; j<Ny_cell; ++j)
         for(int k=0; k<Nz_cell; ++k)
      {
         conserved[i][j][k].density  = primitive[i+ng_cell][j+ng_cell][k+ng_cell].density;
         conserved[i][j][k].momentum.equ(primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity,primitive[i+ng_cell][j+ng_cell][k+ng_cell].density);
         conserved[i][j][k].energy = 0.5*primitive[i+ng_cell][j+ng_cell][k+ng_cell].density*primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity.square() 
                           + primitive[i+ng_cell][j+ng_cell][k+ng_cell].pressure/(param.GAMMA-1.0);
      }
}

//------------------------------------------------------------------------------
// Convert primitive to entropy
//------------------------------------------------------------------------------
void FV::prim_to_ent ()
{

   for(int i=ng_cell; i<Nx_cell+ng_cell; ++i)
      for(int j=ng_cell; j<Ny_cell+ng_cell; ++j)
         for(int k=ng_cell; k<Nz_cell+ng_cell; ++k)
		 {
			double beta = primitive[i][j][k].density/(2*primitive[i][j][k].pressure);
			double s = Entropy (primitive[i][j][k]);
			entropy_var[i][j][k].entf  = (param.GAMMA - s)/(param.GAMMA - 1) - (primitive[i][j][k].velocity.square())*beta;
			entropy_var[i][j][k].entv.equ(primitive[i][j][k].velocity,(2*beta));
			entropy_var[i][j][k].entl = -2*beta;
		 }
   
   // X-Left ghost cell layer
   for(int i=0; i<ng_cell; ++i)
      for(int j=ng_cell; j<Ny_cell+ng_cell; ++j)
         for(int k=ng_cell; k<Nz_cell+ng_cell; ++k)
		 {
			double beta = primitive[i][j][k].density/(2*primitive[i][j][k].pressure);
			double s = Entropy (primitive[i][j][k]);
			entropy_var[i][j][k].entf  = (param.GAMMA - s)/(param.GAMMA - 1) - (primitive[i][j][k].velocity.square())*beta;
			entropy_var[i][j][k].entv.equ(primitive[i][j][k].velocity,(2*beta));
			entropy_var[i][j][k].entl = -2*beta;
		 }
   
   // Right ghost cell layer
   for(int i=Nx_cell+ng_cell; i<Nx_cell+2*ng_cell; ++i)
      for(int j=ng_cell; j<Ny_cell+ng_cell; ++j)
         for(int k=ng_cell; k<Nz_cell+ng_cell; ++k)
		 {
			double beta = primitive[i][j][k].density/(2*primitive[i][j][k].pressure);
			double s = Entropy (primitive[i][j][k]);
			entropy_var[i][j][k].entf  = (param.GAMMA - s)/(param.GAMMA - 1) - (primitive[i][j][k].velocity.square())*beta;
			entropy_var[i][j][k].entv.equ(primitive[i][j][k].velocity,(2*beta));
			entropy_var[i][j][k].entl = -2*beta;
		 }
   
   // Y-left ghost cell layer
   for(int i=ng_cell; i<Nx_cell+ng_cell; ++i)
      for(int j=0; j<ng_cell; ++j)
         for(int k=ng_cell; k<Nz_cell+ng_cell; ++k)
		 {
			double beta = primitive[i][j][k].density/(2*primitive[i][j][k].pressure);
			double s = Entropy (primitive[i][j][k]);
			entropy_var[i][j][k].entf  = (param.GAMMA - s)/(param.GAMMA - 1) - (primitive[i][j][k].velocity.square())*beta;
			entropy_var[i][j][k].entv.equ(primitive[i][j][k].velocity,(2*beta));
			entropy_var[i][j][k].entl = -2*beta;
		 }
   
   // Y-Right ghost cell layer
   for(int i=ng_cell; i<Nx_cell+ng_cell; ++i)
      for(int j=Ny_cell+ng_cell; j<Ny_cell+2*ng_cell; ++j)
         for(int k=ng_cell; k<Nz_cell+ng_cell; ++k)
		 {
			double beta = primitive[i][j][k].density/(2*primitive[i][j][k].pressure);
			double s = Entropy (primitive[i][j][k]);
			entropy_var[i][j][k].entf  = (param.GAMMA - s)/(param.GAMMA - 1) - (primitive[i][j][k].velocity.square())*beta;
			entropy_var[i][j][k].entv.equ(primitive[i][j][k].velocity,(2*beta));
			entropy_var[i][j][k].entl = -2*beta;
		 
		 }
		 
   // Z-left ghost cell layer
   for(int i=ng_cell; i<Nx_cell+ng_cell; ++i)
      for(int j=ng_cell; j<Ny_cell+ng_cell; ++j)
         for(int k=0; k<ng_cell; ++k)
		 {
			double beta = primitive[i][j][k].density/(2*primitive[i][j][k].pressure);
			double s = Entropy (primitive[i][j][k]);
			entropy_var[i][j][k].entf  = (param.GAMMA - s)/(param.GAMMA - 1) - (primitive[i][j][k].velocity.square())*beta;
			entropy_var[i][j][k].entv.equ(primitive[i][j][k].velocity,(2*beta));
			entropy_var[i][j][k].entl = -2*beta;
		 }
   
   // Z-Right ghost cell layer
   for(int i=ng_cell; i<Nx_cell+ng_cell; ++i)
      for(int j=ng_cell; j<Ny_cell+ng_cell; ++j)
         for(int k=Nz_cell+ng_cell; k<Nz_cell+2*ng_cell; ++k)   
		 {
			double beta = primitive[i][j][k].density/(2*primitive[i][j][k].pressure);
			double s = Entropy (primitive[i][j][k]);
			entropy_var[i][j][k].entf  = (param.GAMMA - s)/(param.GAMMA - 1) - (primitive[i][j][k].velocity.square())*beta;
			entropy_var[i][j][k].entv.equ(primitive[i][j][k].velocity,(2*beta));
			entropy_var[i][j][k].entl = -2*beta;
		 
		 }		 
   
   
    // Wedging edges also needed for viscous flux
    // For the time being, it is set assumin all bc are periodic Need to fix this.
    for(int i=0;i<Nx_cell; ++i)
    {
         entropy_var[ng_cell + i][ng_cell - 1][ng_cell - 1]             = entropy_var[ng_cell + i][ng_cell + Ny_cell - 1][ng_cell + Nz_cell - 1];
         entropy_var[ng_cell + i][ng_cell + Ny_cell][ng_cell- 1]        = entropy_var[ng_cell + i][ng_cell][ng_cell + Nz_cell - 1];
         entropy_var[ng_cell + i][ng_cell - 1][ng_cell + Nz_cell]       = entropy_var[ng_cell + i][ng_cell + Ny_cell - 1][ng_cell];
         entropy_var[ng_cell + i][ng_cell + Ny_cell][ng_cell + Nz_cell] = entropy_var[ng_cell + i][ng_cell][ng_cell];
    }
    for(int j=0;j<Ny_cell; ++j)
    {
         entropy_var[ng_cell - 1][ng_cell + j][ng_cell - 1]             = entropy_var[ng_cell + Nx_cell - 1][ng_cell + j][ng_cell + Nz_cell - 1];
         entropy_var[ng_cell + Nx_cell][ng_cell + j][ng_cell - 1]       = entropy_var[ng_cell][ng_cell + j][ng_cell + Nz_cell - 1];
         entropy_var[ng_cell - 1][ng_cell + j][ng_cell + Nz_cell]       = entropy_var[ng_cell + Nx_cell - 1][ng_cell + j][ng_cell];
         entropy_var[ng_cell + Nx_cell][ng_cell + j][ng_cell + Nz_cell] = entropy_var[ng_cell][ng_cell + j][ng_cell];
    }     
    for(int k=0;k<Nz_cell; ++k)
    {
         entropy_var[ng_cell - 1][ng_cell - 1][ng_cell + k]             = entropy_var[ng_cell + Nx_cell - 1][ng_cell + Ny_cell - 1][ng_cell + k];
         entropy_var[ng_cell + Nx_cell][ng_cell - 1][ng_cell + k]       = entropy_var[ng_cell][ng_cell + Ny_cell - 1][ng_cell + k];
         entropy_var[ng_cell - 1][ng_cell + Ny_cell][ng_cell + k]       = entropy_var[ng_cell + Nx_cell - 1][ng_cell][ng_cell + k];
         entropy_var[ng_cell + Nx_cell][ng_cell + Ny_cell][ng_cell + k] = entropy_var[ng_cell][ng_cell][ng_cell + k];
    }   
   
   
   // Corner also needed for viscous flux. For the time being, it is set assumin al bc are periodic Need to fix this.
    // -1,-1,-1
	entropy_var[ng_cell -1][ng_cell -1][ng_cell -1] = entropy_var[ng_cell + Nx_cell -1][ng_cell + Ny_cell -1][ng_cell + Nz_cell -1];
    
    entropy_var[ng_cell + Nx_cell][ng_cell -1][ng_cell -1] = entropy_var[ng_cell][ng_cell + Ny_cell -1][ng_cell + Nz_cell -1];
    
    entropy_var[ng_cell -1][ng_cell + Ny_cell][ng_cell -1] = entropy_var[ng_cell + Nx_cell -1][ng_cell][ng_cell + Nz_cell -1];
    
    entropy_var[ng_cell -1][ng_cell -1][ng_cell + Nz_cell] = entropy_var[ng_cell + Nx_cell -1][ng_cell + Ny_cell -1][ng_cell];
    
    entropy_var[ng_cell + Nx_cell][ng_cell + Ny_cell][ng_cell -1] = entropy_var[ng_cell][ng_cell][ng_cell + Nz_cell -1];
    
    entropy_var[ng_cell + Nx_cell][ng_cell -1][ng_cell + Nz_cell] = entropy_var[ng_cell][ng_cell + Ny_cell -1][ng_cell];
    
    entropy_var[ng_cell -1][ng_cell + Ny_cell][ng_cell + Nz_cell] = entropy_var[ng_cell + Nx_cell -1][ng_cell][ng_cell];
    
    entropy_var[ng_cell + Nx_cell][ng_cell + Ny_cell][ng_cell + Nz_cell] = entropy_var[ng_cell][ng_cell][ng_cell];
}

//------------------------------------------------------------------------------
// Physical entropy
//------------------------------------------------------------------------------
double FV::Entropy (const PrimVar& state) const
{
   return log(state.pressure) - param.GAMMA*log(state.density);
}

//------------------------------------------------------------------------------
// Enthalpy
//------------------------------------------------------------------------------
double FV::Enthalpy (const PrimVar& state) const
{
   double a = sqrt(param.GAMMA*state.pressure/state.density);
   return a*a/(param.GAMMA-1.0) + 0.5*state.velocity.square();
}

// ------------------------------------------------------------------------------
// Compute tau and q ( at face centers generally)
// ------------------------------------------------------------------------------
// void FV::compute_tau_q (vector<double>& left,
//                                 vector<double>& right)
// {
//   double T_left  = temperature(left);
//   double T_right = temperature(right);
//   double T   = 0.5 * (T_left + T_right);
//   double mu  = viscosity (T);
//   double k   = mu * Cp / param.Pr;
// 
//   tau = (4.0 * mu / 3.0 ) * (right[1] - left[1]) / dx;
//   q= -k * (T_right - T_left) / dx;
//   
// }

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FV::compute_residual (const int rk)
{
   for(int i=0; i<Nx_cell; ++i)
   {   for(int j=0; j<Ny_cell; ++j)
          for(int k=0; k<Nz_cell; ++k) 
             residual[i][j][k].zero();
   }
   
   compute_inv_residual(rk);
   
   if(param.model == ns_es || param.model == ns_kep)
      compute_visc_residual(rk);   
}

//------------------------------------------------------------------------------
// Compute inviscid residual
//------------------------------------------------------------------------------
void FV::compute_inv_residual (const int rk)
{
   for(int i=0; i<Nx_cell; ++i)
   {   for(int j=0; j<Ny_cell; ++j)
          for(int k=0; k<Nz_cell; ++k) 
             residual[i][j][k].zero();
   }
   
   Flux flux;
   vector<EntVar> V(2*ng_cell);
   vector<PrimVar> P(2*ng_cell);
   Vector unit_normal;
   
   // x flux
   unit_normal = 0;
   unit_normal.x = 1;
   // Through X-Left boundary
   

   for(int j=0;j<Ny_cell; ++j)
      for(int k=0; k<Nz_cell; ++k)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			V[l] = entropy_var[l][j+ng_cell][k+ng_cell];
			P[l] = primitive[l][j+ng_cell][k+ng_cell];
		 }
		 num_flux(P,V,unit_normal,flux);
		 residual[0][j][k] -= flux;
	  }
   
   // Interior x-flux
   for(int i=1;i<Nx_cell; ++i)
      for(int j=0;j<Ny_cell; ++j)
         for(int k=0; k<Nz_cell; ++k)
		 {
			for(int l=0; l<2*ng_cell; ++l)
			{
			   V[l] = entropy_var[i+l][j+ng_cell][k+ng_cell];
			   P[l] = primitive[i+l][j+ng_cell][k+ng_cell];
			}
			num_flux(P,V,unit_normal,flux);
			residual[i-1][j][k] += flux;
			residual[i][j][k]   -= flux;
		 }
   
   // X-Right boundary
   for(int j=0;j<Ny_cell; ++j)
	  for(int k=0;k<Nz_cell; ++k)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			V[l] = entropy_var[Nx_cell+l][j+ng_cell][k+ng_cell];
			P[l] = primitive[Nx_cell+l][j+ng_cell][k+ng_cell];
		 }
		 num_flux(P,V,unit_normal,flux);
		 residual[Nx_cell-1][j][k] += flux;
	  }
   
   // y flux: 
   unit_normal = 0;
   unit_normal.y = 1;
   // Y-Left boundary
   for(int i=0;i<Nx_cell; ++i)
      for(int k=0;k<Nz_cell; ++k)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			V[l] = entropy_var[i+ng_cell][l][k+ng_cell];
			P[l] = primitive[i+ng_cell][l][k+ng_cell];
		 }
		 num_flux(P,V,unit_normal,flux);
		 residual[i][0][k] -= flux;
	  }
   
   // Interior Y-flux
   for(int i=0;i<Nx_cell; ++i)
      for(int j=1;j<Ny_cell; ++j)
		 for(int k=0;k<Nz_cell; ++k)   
		 {
			for(int l=0; l<2*ng_cell; ++l)
			{
			   V[l] = entropy_var[i+ng_cell][j+l][k+ng_cell];
			   P[l] = primitive[i+ng_cell][j+l][k+ng_cell];
			}
			num_flux(P,V,unit_normal,flux);
			residual[i][j-1][k] += flux;
			residual[i][j][k]   -= flux;
		 }
   
   // Y-Right boundary
   for(int i=0;i<Nx_cell; ++i)
	  for(int k=0;k<Nz_cell; ++k)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			V[l] = entropy_var[i+ng_cell][Ny_cell+l][k+ng_cell];
			P[l] = primitive[i+ng_cell][Ny_cell+l][k+ng_cell];
		 }
		 num_flux(P,V,unit_normal,flux);
		 residual[i][Ny_cell-1][k] += flux;
	  }
   
   // z flux: 
   unit_normal = 0;
   unit_normal.z = 1;
   // Z-Left boundary
   for(int i=0;i<Nx_cell; ++i)
      for(int j=0;j<Ny_cell; ++j)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			V[l] = entropy_var[i+ng_cell][j+ng_cell][l];
			P[l] = primitive[i+ng_cell][j+ng_cell][l];
		 }
		 num_flux(P,V,unit_normal,flux);
		 residual[i][j][0] -= flux;
	  }
   
   // Interior Z-flux
   for(int i=0;i<Nx_cell; ++i)
      for(int j=0;j<Ny_cell; ++j)
		 for(int k=1;k<Nz_cell; ++k)   
		 {
			for(int l=0; l<2*ng_cell; ++l)
			{
			   V[l] = entropy_var[i+ng_cell][j+ng_cell][k+l];
			   P[l] = primitive[i+ng_cell][j+ng_cell][k+l];
			}
			num_flux(P,V,unit_normal,flux);
			residual[i][j][k-1] += flux;
			residual[i][j][k]   -= flux;
		 }
   
   // Z-Right boundary
   for(int i=0;i<Nx_cell; ++i)
	  for(int j=0;j<Ny_cell; ++j)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			V[l] = entropy_var[i+ng_cell][j+ng_cell][Nz_cell+l];
			P[l] = primitive[i+ng_cell][j+ng_cell][Nz_cell+l];
		 }
		 num_flux(P,V,unit_normal,flux);
		 residual[i][j][Nz_cell-1] += flux;
	  }
   
   // Doing this for ever RK step but using it only for rk==0
   KE_diss_pr = 0.0; // At cell-centers
   for(int i=0; i<Nx_cell; ++i)
      for(int j=0; j<Ny_cell; ++j)
          for(int k=0; k<Nz_cell; ++k) 
          {
             // negative sign for dissipation rate
             KE_diss_pr -= 0.5*primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity.square()*residual[i][j][k].mass_flux
                           +primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity*residual[i][j][k].momentum_flux;
          }
       
    
}


//------------------------------------------------------------------------------
// Choose viscous flux type
//------------------------------------------------------------------------------
void FV::nodal_vlux (int i, int j, int k,
                     Flux &xflux, Flux &yflux, Flux &zflux) const
{
   if(param.model == ns_es)
      nodal_vlux_ES(i,j,k,xflux,yflux,zflux);
   else
      nodal_vlux_KE(i,j,k,xflux,yflux,zflux);   
}


//------------------------------------------------------------------------------
// Compute viscous residual
//------------------------------------------------------------------------------
void FV::compute_visc_residual (const int rk)
{
   Flux xflux, yflux, zflux;
   PrimVar prim_avg;
   Vector Txyz_vel, Tyxz_vel, Tzxy_vel;
   double wt;
   KE_diss_stress = 0.0;

   // Loop over interior nodes of active mesh
   // The viscous flux is subtracted not added. The sign has alread been 
   // acdomadated for in the flux evaluation
   // (i,j,k) here correpsonds to interface (i+h/2,j+h/2,k+h/2)
   for(int i=0; i<Nx_cell-1; ++i)
   {
      for(int j=0; j<Ny_cell-1; ++j)
      {
         for(int k=0; k<Nz_cell-1; ++k)
         {
            nodal_vlux(i,j,k,xflux,yflux,zflux);
            for(int l = 0; l<2;++l)
                for(int m =0; m<2; ++m)
                    for(int n=0; n<2; ++n)
                        residual[i+l][j+m][k+n].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,pow(-1,m)/4.0,pow(-1,n)/4.0);
                
                // Doing this for every RK step but saving it only for the first
                wt = node_w[i+1][j+1][k+1];// Node the index shift in find_node_weight.cc
                nodal_avg_and_vel_derivatives (i, j, k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
                // sign for dissipation alrady include in flux
                KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;
                

		 }   
      }
   } 

   // X-Right boundary
   for(int j=0;j<Ny_cell-1; ++j)
	  for(int k=0;k<Nz_cell-1; ++k)
	  {
         nodal_vlux(Nx_cell-1,j,k,xflux,yflux,zflux);
         for(int m =0; m<2; ++m)
            for(int n=0; n<2; ++n)
                residual[Nx_cell-1][j+m][k+n].sadd(xflux,yflux,zflux,1.0/4.0,pow(-1,m)/4.0,pow(-1,n)/4.0);
                
         wt = node_w[Nx_cell][j+1][k+1];// Node the index shift in find_node_weight.cc
         nodal_avg_and_vel_derivatives (Nx_cell - 1, j, k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
         KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;
	  }
   for(int j=0;j<Ny_cell-1; ++j)
   {
      nodal_vlux(Nx_cell-1,j,-1,xflux,yflux,zflux);
      for(int m =0; m<2; ++m)
         residual[Nx_cell-1][j+m][0].sadd(xflux,yflux,zflux,1.0/4.0,pow(-1,m)/4.0,-1.0/4.0);
      
      wt = node_w[Nx_cell][j+1][0];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (Nx_cell - 1, j, -1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;      
      
      nodal_vlux(Nx_cell-1,j,Nz_cell-1,xflux,yflux,zflux);
      for(int m =0; m<2; ++m)
         residual[Nx_cell-1][j+m][Nz_cell-1].sadd(xflux,yflux,zflux,1.0/4.0,pow(-1,m)/4.0,1.0/4.0);  
      wt = node_w[Nx_cell][j+1][Nz_cell];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (Nx_cell - 1, j, Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;       
   } 
   for(int k=0;k<Nz_cell-1; ++k)
   {
      nodal_vlux(Nx_cell-1,-1,k,xflux,yflux,zflux);
      for(int n =0; n<2; ++n)
         residual[Nx_cell-1][0][k+n].sadd(xflux,yflux,zflux,1.0/4.0,-1.0/4.0,pow(-1,n)/4.0);
      wt = node_w[Nx_cell][0][k+1];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (Nx_cell-1,-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;      
      
      nodal_vlux(Nx_cell-1,Ny_cell-1,k,xflux,yflux,zflux);
      for(int n =0; n<2; ++n)
         residual[Nx_cell-1][Ny_cell-1][k+n].sadd(xflux,yflux,zflux,1.0/4.0,1.0/4.0,pow(-1,n)/4.0);
      wt = node_w[Nx_cell][Ny_cell][k+1];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (Nx_cell-1,Ny_cell-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;    
   }     
       
      
      
   // X-Left boundary
   for(int j=0;j<Ny_cell-1; ++j)
	  for(int k=0;k<Nz_cell-1; ++k)
	  {
         nodal_vlux(-1,j,k,xflux,yflux,zflux);
         for(int m =0; m<2; ++m)
            for(int n=0; n<2; ++n)
                residual[0][j+m][k+n].sadd(xflux,yflux,zflux,-1.0/4.0,pow(-1,m)/4.0,pow(-1,n)/4.0);
         wt = node_w[0][j+1][k+1];// Node the index shift in find_node_weight.cc
         nodal_avg_and_vel_derivatives (- 1, j, k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
         KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;       
	  } 
   for(int j=0;j<Ny_cell-1; ++j)
   {
      nodal_vlux(-1,j,-1,xflux,yflux,zflux);
      for(int m =0; m<2; ++m)
         residual[0][j+m][0].sadd(xflux,yflux,zflux,-1.0/4.0,pow(-1,m)/4.0,-1.0/4.0);
      wt = node_w[0][j+1][0];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (-1,j,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;       
      
      nodal_vlux(-1,j,Nz_cell-1,xflux,yflux,zflux);
      for(int m =0; m<2; ++m)
         residual[0][j+m][Nz_cell-1].sadd(xflux,yflux,zflux,-1.0/4.0,pow(-1,m)/4.0,1.0/4.0);   
      wt = node_w[0][j+1][Nz_cell];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (-1,j,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;       
   } 
   for(int k=0;k<Nz_cell-1; ++k)
   {
      nodal_vlux(-1,-1,k,xflux,yflux,zflux);
      for(int n =0; n<2; ++n)
         residual[0][0][k+n].sadd(xflux,yflux,zflux,-1.0/4.0,-1.0/4.0,pow(-1,n)/4.0);
      wt = node_w[0][0][k+1];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (-1,-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;     
      
      nodal_vlux(-1,Ny_cell-1,k,xflux,yflux,zflux);
      for(int n =0; n<2; ++n)
         residual[0][Ny_cell-1][k+n].sadd(xflux,yflux,zflux,-1.0/4.0,1.0/4.0,pow(-1,n)/4.0);
      wt = node_w[0][Ny_cell][k+1];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (-1,Ny_cell-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;     
   }     
      
      
   // Y-Right boundary
   for(int i=0;i<Nx_cell-1; ++i)
	  for(int k=0;k<Nz_cell-1; ++k)
	  {
         nodal_vlux(i,Ny_cell-1,k,xflux,yflux,zflux);
         for(int l = 0; l<2;++l)
            for(int n=0; n<2; ++n)
                residual[i+l][Ny_cell -1][k+n].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,1.0/4.0,pow(-1,n)/4.0);
         wt = node_w[i+1][Ny_cell][k+1];// Node the index shift in find_node_weight.cc
         nodal_avg_and_vel_derivatives (i,Ny_cell-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
         KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;       
	  }
   for(int i=0;i<Nx_cell-1; ++i)
   {
      nodal_vlux(i,Ny_cell-1,-1,xflux,yflux,zflux);
      for(int l = 0; l<2;++l)
         residual[i+l][Ny_cell -1][0].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,1.0/4.0,-1.0/4.0);
      wt = node_w[i+1][Ny_cell][0];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (i,Ny_cell-1,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;    
         
      nodal_vlux(i,Ny_cell-1,Nz_cell-1,xflux,yflux,zflux);
      for(int l = 0; l<2;++l)
         residual[i+l][Ny_cell -1][Nz_cell-1].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,1.0/4.0,1.0/4.0);   
      wt = node_w[i+1][Ny_cell][Nz_cell];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (i,Ny_cell-1,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;       
   }  
   for(int k=0;k<Nz_cell-1; ++k)
   {
      nodal_vlux(-1,Ny_cell-1,k,xflux,yflux,zflux);
      for(int n = 0; n<2;++n)
         residual[0][Ny_cell -1][k+n].sadd(xflux,yflux,zflux,-1.0/4.0,1.0/4.0,pow(-1,n)/4.0);
      wt = node_w[0][Ny_cell][k+1];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (-1,Ny_cell-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;          
         
      nodal_vlux(Nx_cell-1,Ny_cell-1,k,xflux,yflux,zflux);
      for(int n = 0; n<2;++n)
         residual[Nx_cell-1][Ny_cell -1][k+n].sadd(xflux,yflux,zflux,1.0/4.0,1.0/4.0,pow(-1,n)/4.0);
      wt = node_w[Nx_cell][Ny_cell][k+1];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (Nx_cell-1,Ny_cell-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;     
   } 
      
   // Y-Left boundary
   for(int i=0;i<Nx_cell-1; ++i)
	  for(int k=0;k<Nz_cell-1; ++k)
	  {
         nodal_vlux(i,-1,k,xflux,yflux,zflux);
         for(int l = 0; l<2;++l)
            for(int n=0; n<2; ++n)
                residual[i+l][0][k+n].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,-1.0/4.0,pow(-1,n)/4.0);
         wt = node_w[i+1][0][k+1];// Node the index shift in find_node_weight.cc
         nodal_avg_and_vel_derivatives (i,-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
         KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;             
	  }   
   for(int i=0;i<Nx_cell-1; ++i)
   {
      nodal_vlux(i,-1,-1,xflux,yflux,zflux);
      for(int l = 0; l<2;++l)
         residual[i+l][0][0].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,-1.0/4.0,-1.0/4.0);
      wt = node_w[i+1][0][0];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (i,-1,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;      
         
      nodal_vlux(i,-1,Nz_cell-1,xflux,yflux,zflux);
      for(int l = 0; l<2;++l)
         residual[i+l][0][Nz_cell-1].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,-1.0/4.0,1.0/4.0);   
      wt = node_w[i+1][0][Nz_cell];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (i,-1,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;     
   }  
   for(int k=0;k<Nz_cell-1; ++k)
   {
      nodal_vlux(-1,-1,k,xflux,yflux,zflux);
      for(int n = 0; n<2;++n)
         residual[0][0][k+n].sadd(xflux,yflux,zflux,-1.0/4.0,-1.0/4.0,pow(-1,n)/4.0);
      wt = node_w[0][0][k+1];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (-1,-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;     
         
      nodal_vlux(Nx_cell-1,-1,k,xflux,yflux,zflux);
      for(int n = 0; n<2;++n)
         residual[Nx_cell-1][0][k+n].sadd(xflux,yflux,zflux,1.0/4.0,-1.0/4.0,pow(-1,n)/4.0);
      wt = node_w[Nx_cell][0][k+1];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (Nx_cell-1,-1,k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;    
   }     
      
   // Z-Right boundary
   for(int i=0;i<Nx_cell-1; ++i)
	  for(int j=0;j<Ny_cell-1; ++j)
	  {
         nodal_vlux(i,j,Nz_cell-1,xflux,yflux,zflux);
         for(int l = 0; l<2;++l)
            for(int m=0; m<2; ++m)
                residual[i+l][j+m][Nz_cell -1].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,pow(-1,m)/4.0,1.0/4.0);
         wt = node_w[i+1][j+1][Nz_cell];// Node the index shift in find_node_weight.cc
         nodal_avg_and_vel_derivatives (i,j,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
         KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;              
	  }
      
   for(int i=0;i<Nx_cell-1; ++i)
   {
      nodal_vlux(i,-1,Nz_cell-1,xflux,yflux,zflux);
      for(int l = 0; l<2;++l)
         residual[i+l][0][Nz_cell -1].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,-1.0/4.0,1.0/4.0);
      wt = node_w[i+1][0][Nz_cell];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (i,-1,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;  
         
      nodal_vlux(i,Ny_cell-1,Nz_cell-1,xflux,yflux,zflux);
      for(int l = 0; l<2;++l)
         residual[i+l][Ny_cell-1][Nz_cell -1].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,1.0/4.0,1.0/4.0);   
      wt = node_w[i+1][Ny_cell][Nz_cell];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (i,Ny_cell-1,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;     
   }  
   
   for(int j=0;j<Ny_cell-1; ++j)
   {
      nodal_vlux(-1,j,Nz_cell-1,xflux,yflux,zflux);
         for(int m=0; m<2; ++m)
            residual[0][j+m][Nz_cell -1].sadd(xflux,yflux,zflux,-1.0/4.0,pow(-1,m)/4.0,1.0/4.0);
      wt = node_w[0][j+1][Nz_cell];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (-1,j,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;         
            
      nodal_vlux(Nx_cell-1,j,Nz_cell-1,xflux,yflux,zflux);
         for(int m=0; m<2; ++m)
            residual[Nx_cell-1][j+m][Nz_cell -1].sadd(xflux,yflux,zflux,1.0/4.0,pow(-1,m)/4.0,1.0/4.0);  
      wt = node_w[Nx_cell][j+1][Nz_cell];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (Nx_cell-1,j,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;       
               
   } 
      
   // Z-Left boundary
   for(int i=0;i<Nx_cell-1; ++i)
	  for(int j=0;j<Ny_cell-1; ++j)
	  {
         nodal_vlux(i,j,-1,xflux,yflux,zflux);
         for(int l = 0; l<2;++l)
            for(int m=0; m<2; ++m)
                residual[i+l][j+m][0].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,pow(-1,m)/4.0,-1.0/4.0);
         wt = node_w[i+1][j+1][0];// Node the index shift in find_node_weight.cc
         nodal_avg_and_vel_derivatives (i,j,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
         KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;           
	  } 
      
   for(int i=0;i<Nx_cell-1; ++i)
   {
      nodal_vlux(i,-1,-1,xflux,yflux,zflux);
      for(int l = 0; l<2;++l)
         residual[i+l][0][0].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,-1.0/4.0,-1.0/4.0);
      wt = node_w[i+1][0][0];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (i,-1,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;   
         
      nodal_vlux(i,Ny_cell-1,-1,xflux,yflux,zflux);
      for(int l = 0; l<2;++l)
         residual[i+l][Ny_cell-1][0].sadd(xflux,yflux,zflux,pow(-1,l)/4.0,1.0/4.0,-1.0/4.0);   
      wt = node_w[i+1][Ny_cell][0];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (i,Ny_cell-1,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;    
   }  
   
   for(int j=0;j<Ny_cell-1; ++j)
   {
      nodal_vlux(-1,j,-1,xflux,yflux,zflux);
         for(int m=0; m<2; ++m)
            residual[0][j+m][0].sadd(xflux,yflux,zflux,-1.0/4.0,pow(-1,m)/4.0,-1.0/4.0);
      wt = node_w[0][j+1][0];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (-1,j,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;          
            
      nodal_vlux(Nx_cell-1,j,-1,xflux,yflux,zflux);
         for(int m=0; m<2; ++m)
            residual[Nx_cell-1][j+m][0].sadd(xflux,yflux,zflux,1.0/4.0,pow(-1,m)/4.0,-1.0/4.0);      
      wt = node_w[Nx_cell][j+1][0];// Node the index shift in find_node_weight.cc
      nodal_avg_and_vel_derivatives (Nx_cell-1,j,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
      KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;         
   }   
   
   //Corner fluxes (8 of them)
   // (0,0,0)
   nodal_vlux(-1,-1,-1,xflux,yflux,zflux);
   residual[0][0][0].sadd(xflux,yflux,zflux,-1.0/4.0,-1.0/4.0,-1.0/4.0);
   wt = node_w[0][0][0];// Node the index shift in find_node_weight.cc
   nodal_avg_and_vel_derivatives (-1,-1,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;   
   
   // (N,0,0)
   nodal_vlux(Nx_cell-1,-1,-1,xflux,yflux,zflux);
   residual[Nx_cell-1][0][0].sadd(xflux,yflux,zflux,1.0/4.0,-1.0/4.0,-1.0/4.0);
   wt = node_w[Nx_cell][0][0];// Node the index shift in find_node_weight.cc
   nodal_avg_and_vel_derivatives (Nx_cell-1,-1,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;   
   
   // (0,N,0)
   nodal_vlux(-1,Ny_cell-1,-1,xflux,yflux,zflux);
   residual[0][Ny_cell-1][0].sadd(xflux,yflux,zflux,-1.0/4.0,1.0/4.0,-1.0/4.0);
   wt = node_w[0][Ny_cell][0];// Node the index shift in find_node_weight.cc
   nodal_avg_and_vel_derivatives (-1,Ny_cell-1,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;  
   
   // (0,0,N)
   nodal_vlux(-1,-1,Nz_cell-1,xflux,yflux,zflux);
   residual[0][0][Nz_cell-1].sadd(xflux,yflux,zflux,-1.0/4.0,-1.0/4.0,1.0/4.0);
   wt = node_w[0][0][Nz_cell];// Node the index shift in find_node_weight.cc
   nodal_avg_and_vel_derivatives (-1,-1,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;  
   
   // (N,N,0)
   nodal_vlux(Nx_cell-1,Ny_cell-1,-1,xflux,yflux,zflux);
   residual[Nx_cell-1][Ny_cell-1][0].sadd(xflux,yflux,zflux,1.0/4.0,1.0/4.0,-1.0/4.0);
   wt = node_w[Nx_cell][Ny_cell][0];// Node the index shift in find_node_weight.cc
   nodal_avg_and_vel_derivatives (Nx_cell-1,Ny_cell-1,-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;   
   
   // (N,0,N)
   nodal_vlux(Nx_cell-1,-1,Nz_cell-1,xflux,yflux,zflux);
   residual[Nx_cell-1][0][Nz_cell-1].sadd(xflux,yflux,zflux,1.0/4.0,-1.0/4.0,1.0/4.0);
   wt = node_w[Nx_cell][0][Nz_cell];// Node the index shift in find_node_weight.cc
   nodal_avg_and_vel_derivatives (Nx_cell-1,-1,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;   
   
   // (0,N,N)
   nodal_vlux(-1,Ny_cell-1,Nz_cell-1,xflux,yflux,zflux);
   residual[0][Ny_cell-1][Nz_cell-1].sadd(xflux,yflux,zflux,-1.0/4.0,1.0/4.0,1.0/4.0);
   wt = node_w[0][Nz_cell][Nz_cell];// Node the index shift in find_node_weight.cc
   nodal_avg_and_vel_derivatives (-1,Nz_cell-1,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;  
   
   // (N,N,N)
   nodal_vlux(Nx_cell-1,Ny_cell-1,Nz_cell-1,xflux,yflux,zflux);
   residual[Nx_cell-1][Ny_cell-1][Nz_cell-1].sadd(xflux,yflux,zflux,1.0/4.0,1.0/4.0,1.0/4.0);
   wt = node_w[Nx_cell][Nz_cell][Nz_cell];// Node the index shift in find_node_weight.cc
   nodal_avg_and_vel_derivatives (Nx_cell-1,Nz_cell-1,Nz_cell-1, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   KE_diss_stress -= (Txyz_vel*xflux.momentum_flux*dx + Tyxz_vel*yflux.momentum_flux*dy + Tzxy_vel*zflux.momentum_flux*dz)*wt;  
}


//------------------------------------------------------------------------------
// Update solution
//------------------------------------------------------------------------------
void FV::update_solution (const int rk)
{
   double factor = dt/(dx*dy*dz);
   if (param.time_scheme == rk1 || param.time_scheme == ssprk3)
   {
      for(int i=0; i<Nx_cell; ++i)
         for(int j=0; j<Ny_cell; ++j)
            for(int k=0; k<Nz_cell; ++k)
			   {
               conserved[i][j][k] *= brks[rk];
               conserved[i][j][k].sadd(conserved_old[i][j][k],residual[i][j][k],arks[rk],-factor*brks[rk]);
			   }                     
   }       
   else if (param.time_scheme == jameson_rk4)
   {
      for(int i=0; i<Nx_cell; ++i)
         for(int j=0; j<Ny_cell; ++j)
            for(int k=0; k<Nz_cell; ++k)
			{
			   conserved[i][j][k].equ(conserved_old[i][j][k], residual[i][j][k],1.0,-factor*jameson_rks[rk]); 
			}                       
   }
   else if (param.time_scheme == rk4)
   {
      for(int i=0; i<Nx_cell; ++i)
         for(int j=0; j<Ny_cell; ++j)
            for(int k=0; k<Nz_cell; ++k)
			{
			   residual2[i][j][k].sadd(residual[i][j][k],1.0/rk4_rks[rk]);
			   if(rk <3)
			      conserved[i][j][k].equ(conserved_old[i][j][k], residual[i][j][k],1.0,-factor*rk4_rks[rk+1]); 
			   else
			      conserved[i][j][k].equ(conserved_old[i][j][k], residual2[i][j][k],1.0,-factor/6.0);    
			}                       
   }                           
}


//------------------------------------------------------------------------------
// compute norm of residual
//------------------------------------------------------------------------------
void FV::compute_residual_norm ()
{
   res_norm.zero();
   
   for(int i=0; i<Nx_cell; ++i)
	 for(int j=0; j<Ny_cell; ++j)
	    for(int k=0; k<Nz_cell; ++k)
		{
		   res_norm.mass_flux       += pow(residual[i][j][k].mass_flux/dx/dy/dz,2);
		   res_norm.momentum_flux.x += pow(residual[i][j][k].momentum_flux.x/dx/dy/dz,2);
		   res_norm.momentum_flux.y += pow(residual[i][j][k].momentum_flux.y/dx/dy/dz,2);
		   res_norm.momentum_flux.z += pow(residual[i][j][k].momentum_flux.z/dx/dy/dz,2);
		   res_norm.energy_flux     += pow(residual[i][j][k].energy_flux/dx/dy/dz,2);
		}		                  
   
   res_norm.mass_flux       = sqrt(res_norm.mass_flux);
   res_norm.momentum_flux.x = sqrt(res_norm.momentum_flux.x);
   res_norm.momentum_flux.y = sqrt(res_norm.momentum_flux.y);
   res_norm.momentum_flux.z = sqrt(res_norm.momentum_flux.z);
   res_norm.energy_flux     = sqrt(res_norm.energy_flux);
   
   res_norm*=(dx*dy*dz);
}

//------------------------------------------------------------------------------
// Start the computations
//------------------------------------------------------------------------------
void FV::run ()
{
   for(int rl=0; rl<param.mesh_levels; ++rl)
   {
	  Nx_cell = param.Nx_cell[rl];
	  Ny_cell = param.Ny_cell[rl];
	  Nz_cell = param.Nz_cell[rl];
	  
	  cout<<"\n\n   Solving on mesh with size "<<Nx_cell<<"x"<<Ny_cell<<"x"<<Nz_cell<<endl;
	     
	  initialize();
	  prim_to_con ();
	  prim_to_ent();
      find_vorticity();
	  counter = 0;
	  time = 0.0;
	  iter = 0;
	  
	  double dt_print = param.final_time/100;
	  double t_print_next = dt_print;
	  int t_percent = 0;
      
	  if(print_soln)
		 output();    
         
      double t_interval_size = param.final_time/param.t_stamps;
      double t_soln_save = t_interval_size;   
      
      cout<<"  ---Solving ...\n";
      if(!disp_log)
         cout<<"          "<<t_percent<<"%"<<flush;
	  while (time < param.final_time && iter < param.max_iter)
	  {
		 conserved_old = conserved;
		 compute_dt ();
		 if(time+dt > t_soln_save) dt = t_soln_save - time;
		 
		 if(param.time_scheme==rk4)
	     {
	        for(int i=0; i<Nx_cell; ++i)
			{   for(int j=0; j<Ny_cell; ++j)
				   for(int k=0; k<Nz_cell; ++k) 
					  residual2[i][j][k].zero();
			}
	     }
		 for(int rk=0; rk<param.n_rks; ++rk)
		 {
			prim_to_ent();
			compute_residual (rk);
            if(rk==0)
              compute_global(); // NOT SAVED AT FINAL TIME BECAUSE OF THE WAY DISSIPATION RATE IS EVALUATED
			update_solution (rk);
			con_to_prim ();
		 }
		 time += dt;
		 ++iter;
		 compute_residual_norm();
	  
		 if(disp_log)
			cout << left
				 << "ITER = "<< setw(8) << iter << " "
				 << scientific
				 << setprecision (4)
				 << "TIME = "<<time << " "
				 << "RES_MASS = "<<res_norm.mass_flux << " "
				 << "RES_MOM_X = "<<res_norm.momentum_flux.x << " "
				 << "RES_MOM_y = "<<res_norm.momentum_flux.x << " "
				 << "RES_MOM_z = "<<res_norm.momentum_flux.z << " "
				 << "RES_ENERGY = "<<res_norm.energy_flux << endl;
	     else if(time >= t_print_next)
	     {
	        t_percent++;
	        cout<<"\r          "<<t_percent<<"%"<<flush;		 
	        t_print_next+=dt_print;
	     }
         find_vorticity();
 
		 if(time == t_soln_save)
         {
            if(print_soln) output();
            t_soln_save += t_interval_size;
         }   
	  }
	  if(!disp_log)
         cout<<"\r          "<<++t_percent<<"%\n"<<flush;
	  if(param.OOC)
		 compute_error(rl);
	  clear_data_structs();	 
   }

   if(param.OOC)
      compute_OOC();

}

//------------------------------------------------------------------------------
// Compute temperature given primitive variables
//------------------------------------------------------------------------------
 double FV::Temperature(const PrimVar &prim) const
 {
    return prim.pressure / (param.gas_const * prim.density);
 }

//------------------------------------------------------------------------------
// Compute temperature given primitive variables
//------------------------------------------------------------------------------
// double FV::enthalpy(const vector<double>& prim) const
// {
//    return param.GAMMA * prim[2] / (prim[0] * (param.GAMMA-1.0)) + 
//           0.5 * pow(prim[1], 2);
// }

