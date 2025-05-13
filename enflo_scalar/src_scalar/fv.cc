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
//------------------------------------------------------------------------------
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
   solution.resize (Nx_cell + 2*ng_cell);
   solution_old.resize (Nx_cell + 2*ng_cell);
   residual.resize (Nx_cell);
   if(param.time_scheme == rk4)
      residual2.resize (Nx_cell);
   node_w.resize(Nx_cell + 1);
   
   for(int i=0; i<Nx_cell + 2*ng_cell; ++i)
   {
	  solution[i].resize (Ny_cell + 2*ng_cell, vector<double>(Nz_cell + 2*ng_cell));
	  solution_old[i].resize (Ny_cell + 2*ng_cell, vector<double>(Nz_cell + 2*ng_cell));
	  if(i>=ng_cell && i < Nx_cell + ng_cell)
	  {
		 residual[i-ng_cell].resize (Ny_cell, vector<double>(Nz_cell));
		 if(param.time_scheme == rk4)
            residual2[i-ng_cell].resize (Ny_cell, vector<double>(Nz_cell));
	  }
   }
   
   for(int i=0; i<Nx_cell+1; ++i)
   {
      node_w[i].resize (Ny_cell + 1, vector<double>(Nz_cell+1));
   }  
   find_node_weights();
	  
   for(int i=0; i<Nx_cell; ++i)
   {
      for(int j=0; j<Ny_cell; ++j)
      {
         for(int k=0; k<Nz_cell; ++k)
		 {
			param.initial_condition.value (cell_cc[i][j][k],solution[i+ng_cell][j+ng_cell][k+ng_cell]);
	
		 }   
      }
   }
   
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
   solution.clear();
   residual.clear();
   residual2.clear();
   solution_old.clear();
   
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
			    solution[i][ng_cell + j][ng_cell + k] = solution[Nx_cell + i][ng_cell + j][ng_cell + k];
	}
	else if(param.bc[0] == wall)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				solution[i][ng_cell + j][ng_cell + k] = -solution[2*ng_cell - 1 - i][ng_cell + j][ng_cell + k];
			 }
	}
	else if(param.bc[0] == neumann)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				solution[i][ng_cell + j][ng_cell + k] = solution[ng_cell][ng_cell + j][ng_cell + k];
			 }
	}
	
	// X-Right Wall
	if(param.bc[1] == periodic)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			    solution[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k] = solution[ng_cell + i][ng_cell + j][ng_cell + k];	
	}
	else if(param.bc[1] == wall)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				solution[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k] = - solution[Nx_cell + ng_cell - 1 - i][ng_cell + j][ng_cell + k];
			 }
	}
	else if(param.bc[1] == neumann)
	{
	   for(int i = 0; i< ng_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				solution[Nx_cell + ng_cell + i][ng_cell + j][ng_cell + k] = solution[Nx_cell + ng_cell - 1][ng_cell + j][ng_cell + k];
			 }
	}
    
    // Y-Left Wall
	if(param.bc[2] == periodic)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			    solution[ng_cell + i][j][ng_cell + k] = solution[ng_cell + i][Ny_cell + j][ng_cell + k];
	}
	else if(param.bc[2] == wall)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				solution[ng_cell + i][j][ng_cell + k] = - solution[ng_cell + i][2*ng_cell - 1 - j][ng_cell + k];
			 }
	}
	else if(param.bc[2] == neumann)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				solution[ng_cell + i][j][ng_cell + k] = solution[ng_cell + i][ng_cell][ng_cell + k];
			 }		  
	} 
    
    // Y-Right Wall
	if(param.bc[3] == periodic)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
			 for(int k = 0; k< Nz_cell; ++k)
			    solution[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k] = solution[ng_cell + i][ng_cell + j][ng_cell + k];
	}
	else if(param.bc[3] == wall)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				solution[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k] = - solution[ng_cell + i][Ny_cell + ng_cell - 1 - j][ng_cell + k];
			 }
	}
	else if(param.bc[3] == neumann)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< ng_cell; ++j)
		     for(int k = 0; k< Nz_cell; ++k)
			 {
				solution[ng_cell + i][Ny_cell + ng_cell + j][ng_cell + k] = solution[ng_cell + i][Ny_cell + ng_cell - 1][ng_cell + k];
			 }
		  
	   
	}
	
	
	// Z-Left Wall
	if(param.bc[4] == periodic)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			    solution[ng_cell + i][ng_cell + j][k] = solution[ng_cell + i][ng_cell + j][Nz_cell + k];
	}
	else if(param.bc[4] == wall)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			 {
				solution[ng_cell + i][ng_cell + j][k] = -solution[ng_cell + i][ng_cell + j][2*ng_cell - 1 - k];
			 }
	}
	else if(param.bc[4] == neumann)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			 {
				solution[ng_cell + i][ng_cell + j][k] = solution[ng_cell + i][ng_cell + j][ng_cell];
			 }		  
	} 
	
	// Z-Right Wall
	if(param.bc[5] == periodic)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
			 for(int k = 0; k< ng_cell; ++k)
			    solution[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k] = solution[ng_cell + i][ng_cell + j][ng_cell + k];
	}
	else if(param.bc[5] == wall)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			 {
				solution[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k] = -solution[ng_cell + i][ng_cell + j][Nz_cell + ng_cell - 1 - k];
	         }
	}
	else if(param.bc[5] == neumann)
	{
	   for(int i = 0; i< Nx_cell; ++i)
		  for(int j = 0; j< Ny_cell; ++j)
		     for(int k = 0; k< ng_cell; ++k)
			 {
				solution[ng_cell + i][ng_cell + j][Nz_cell + ng_cell + k] = solution[ng_cell + i][ng_cell + j][Nz_cell + ng_cell - 1];
			 }
	   
	}
    
	
}

//------------------------------------------------------------------------------
// Compute time step (DX and DY factor???)
//------------------------------------------------------------------------------
void FV::compute_dt ()
{
   dt = 1.0e20;
   if(param.model == linadv)
   {
	 for(int i=0; i<Nx_cell; ++i)
		for(int j=0; j<Ny_cell; ++j)
		   for(int k=0; k<Nz_cell; ++k)
		   {
			  double speed = param.lin_vel.norm();
			  dt = min (dt, min(dx,min(dy,dz))/speed);
		   }
   } 
   else
   {
     for(int i=0; i<Nx_cell; ++i)
		for(int j=0; j<Ny_cell; ++j)
		   for(int k=0; k<Nz_cell; ++k)
		   {
			  double speed = fabs(solution[i+ng_cell][j+ng_cell][k+ng_cell]);
			  dt = min (dt, min(dx,min(dy,dz))/speed);
		   }
   }
   dt *= param.cfl;
   
}


//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FV::compute_residual (const int rk)
{
   for(int i=0; i<Nx_cell; ++i)
   {   for(int j=0; j<Ny_cell; ++j)
          for(int k=0; k<Nz_cell; ++k) 
             residual[i][j][k] = 0.0;
   }
   
   double flux;
   vector<double> Soln(2*ng_cell);
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
			Soln[l] = solution[l][j+ng_cell][k+ng_cell];
		 }
		 num_flux(Soln,unit_normal,flux);
		 residual[0][j][k] -= flux;
	  }
   
   // Interior x-flux
   for(int i=1;i<Nx_cell; ++i)
      for(int j=0;j<Ny_cell; ++j)
         for(int k=0; k<Nz_cell; ++k)
		 {
			for(int l=0; l<2*ng_cell; ++l)
			{
			   Soln[l] = solution[i+l][j+ng_cell][k+ng_cell];
			}
			//cout<<cell_cc[i][j].x<<","<<cell_cc[i][j].y<<endl;
			num_flux(Soln,unit_normal,flux);
			residual[i-1][j][k] += flux;
			residual[i][j][k]   -= flux;
		 }
   
   // X-Right boundary
   for(int j=0;j<Ny_cell; ++j)
	  for(int k=0;k<Nz_cell; ++k)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			Soln[l] = solution[Nx_cell+l][j+ng_cell][k+ng_cell];
		 }
		 num_flux(Soln,unit_normal,flux);
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
			Soln[l] = solution[i+ng_cell][l][k+ng_cell];
		 }
		 num_flux(Soln,unit_normal,flux);
		 residual[i][0][k] -= flux;
	  }
   //cout<<"DONE BOTTOM"<<endl;
   
   // Interior Y-flux
   for(int i=0;i<Nx_cell; ++i)
      for(int j=1;j<Ny_cell; ++j)
		 for(int k=0;k<Nz_cell; ++k)   
		 {
			for(int l=0; l<2*ng_cell; ++l)
			{
			   Soln[l] = solution[i+ng_cell][j+l][k+ng_cell];
			}
			num_flux(Soln,unit_normal,flux);
			residual[i][j-1][k] += flux;
			residual[i][j][k]   -= flux;
		 }
   
   // Y-Right boundary
   for(int i=0;i<Nx_cell; ++i)
	  for(int k=0;k<Nz_cell; ++k)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			Soln[l] = solution[i+ng_cell][Ny_cell+l][k+ng_cell];
		 }
		 num_flux(Soln,unit_normal,flux);
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
			Soln[l] = solution[i+ng_cell][j+ng_cell][l];
		 }
		 num_flux(Soln,unit_normal,flux);
		 residual[i][j][0] -= flux;
	  }
   //cout<<"DONE BOTTOM"<<endl;
   
   // Interior Z-flux
   for(int i=0;i<Nx_cell; ++i)
      for(int j=0;j<Ny_cell; ++j)
		 for(int k=1;k<Nz_cell; ++k)   
		 {
			for(int l=0; l<2*ng_cell; ++l)
			{
			   Soln[l] = solution[i+ng_cell][j+ng_cell][k+l];
			}
			num_flux(Soln,unit_normal,flux);
			residual[i][j][k-1] += flux;
			residual[i][j][k]   -= flux;
		 }
   
   // Z-Right boundary
   for(int i=0;i<Nx_cell; ++i)
	  for(int j=0;j<Ny_cell; ++j)
	  {
		 for(int l=0; l<2*ng_cell; ++l)
		 {
			Soln[l] = solution[i+ng_cell][j+ng_cell][Nz_cell+l];
		 }
		 num_flux(Soln,unit_normal,flux);
		 residual[i][j][Nz_cell-1] += flux;
	  }
       
    
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
			   solution[i+ng_cell][j+ng_cell][k+ng_cell] *= brks[rk];
			   solution[i+ng_cell][j+ng_cell][k+ng_cell] +=(solution_old[i+ng_cell][j+ng_cell][k+ng_cell]*arks[rk] -factor*brks[rk]*residual[i][j][k]);
			}                     
   }       
   else if (param.time_scheme == jameson_rk4)
   {
      for(int i=0; i<Nx_cell; ++i)
         for(int j=0; j<Ny_cell; ++j)
            for(int k=0; k<Nz_cell; ++k)
			{
			   solution[i+ng_cell][j+ng_cell][k+ng_cell] = (solution_old[i+ng_cell][j+ng_cell][k+ng_cell] -factor*jameson_rks[rk]*residual[i][j][k]); 
			}                       
   }
   else if (param.time_scheme == rk4)
   {
      for(int i=0; i<Nx_cell; ++i)
         for(int j=0; j<Ny_cell; ++j)
            for(int k=0; k<Nz_cell; ++k)
			{
			   residual2[i][j][k]+=(residual[i][j][k]/rk4_rks[rk]);
			   if(rk <3)
			      solution[i+ng_cell][j+ng_cell][k+ng_cell] = (solution_old[i+ng_cell][j+ng_cell][k+ng_cell] -factor*rk4_rks[rk+1]*residual[i][j][k]); 
			   else
			      solution[i+ng_cell][j+ng_cell][k+ng_cell] = (solution_old[i+ng_cell][j+ng_cell][k+ng_cell] -factor*residual2[i][j][k]/6.0);    
			}                       
   }                           
}


//------------------------------------------------------------------------------
// compute norm of residual
//------------------------------------------------------------------------------
void FV::compute_residual_norm ()
{
   res_norm = 0.0;
   
   for(int i=0; i<Nx_cell; ++i)
	 for(int j=0; j<Ny_cell; ++j)
	    for(int k=0; k<Nz_cell; ++k)
		{
		   res_norm += pow(residual[i][j][k],2);
		}		                  
   
   res_norm = sqrt(res_norm);
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
	  counter = 0;
	  time = 0.0;
	  iter = 0;
	  
	  double dt_print = param.final_time/100;
	  double t_print_next = dt_print;
	  int t_percent = 0;
      
	  if(print_soln)
		 output();    
		 
      compute_global();		 
         
      double t_interval_size = param.final_time/param.t_stamps;
      double t_soln_save = t_interval_size;   
      
      cout<<"  ---Solving ...\n";
      if(!disp_log)
         cout<<"          "<<t_percent<<"%"<<flush;
      
	  while (time < param.final_time && iter < param.max_iter)
	  {
		 solution_old = solution;
		 compute_dt ();
		 if(time+dt > t_soln_save) dt = t_soln_save - time;
		 
		 if(param.time_scheme==rk4)
	     {
	        for(int i=0; i<Nx_cell; ++i)
			{   for(int j=0; j<Ny_cell; ++j)
				   for(int k=0; k<Nz_cell; ++k) 
					  residual2[i][j][k] = 0.0;
			}
	     }
		 for(int rk=0; rk<param.n_rks; ++rk)
		 {
			compute_residual (rk);
			update_solution (rk);
			apply_bc();
		 }
		 time += dt;
		 ++iter;
		 compute_residual_norm();
		 compute_global();
	  
		 if(disp_log)
			cout << left
				 << "ITER = "<< setw(8) << iter << " "
				 << scientific
				 << setprecision (4)
				 << "TIME = "<<time << " "
				 << "RES = "<<res_norm << endl;
	     else if(time >= t_print_next)
	     {
	        t_percent++;
	        cout<<"\r          "<<t_percent<<"%"<<flush;		 
	        t_print_next+=dt_print;
	     }   
 
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


