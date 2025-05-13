#include <cmath>
#include <iomanip>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Computing OOC
//------------------------------------------------------------------------------
void FV::compute_OOC()
{
   PrimVar L1_OOC, L2_OOC, Linf_OOC;
   
   L1_err_file << left << setw(112) << setfill('-') <<""<<endl;
   L1_err_file << left << setw(6) << setfill(' ') << "Nx"
               << left << setw(6) << setfill(' ') << "Ny" 
               << left << setw(6) << setfill(' ') << "Nz" 
			   << left << setw(17)<< setfill(' ') <<"DENSITY_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"
			   << left << setw(17)<< setfill(' ') <<"VELX_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"
			   << left << setw(17)<< setfill(' ') <<"VELY_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"
			   << left << setw(17)<< setfill(' ') <<"VELZ_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"
			   << left << setw(17)<< setfill(' ') <<"PRESSURE_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"<<endl
               << left << setw(112) << setfill('-') << ""<<endl; 
   
   L2_err_file << left << setw(112) << setfill('-') <<""<<endl;
   L2_err_file << left << setw(6) << setfill(' ') << "Nx"
               << left << setw(6) << setfill(' ') << "Ny" 
               << left << setw(6) << setfill(' ') << "Nz"
			   << left << setw(17)<< setfill(' ') <<"DENSITY_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"
			   << left << setw(17)<< setfill(' ') <<"VELX_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"
			   << left << setw(17)<< setfill(' ') <<"VELY_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"
			   << left << setw(17)<< setfill(' ') <<"VELZ_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"
			   << left << setw(17)<< setfill(' ') <<"PRESSURE_ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"<<endl
               << left << setw(112) << setfill('-') << ""<<endl; 
   
   Linf_err_file << left << setw(112) << setfill('-') <<""<<endl;
   Linf_err_file << left << setw(6) << setfill(' ') << "Nx"
				 << left << setw(6) << setfill(' ') << "Ny" 
				 << left << setw(6) << setfill(' ') << "Nz"
				 << left << setw(17)<< setfill(' ') <<"DENSITY_ERR"
				 << left << setw(8) << setfill(' ') <<"OOC"
				 << left << setw(17)<< setfill(' ') <<"VELX_ERR"
				 << left << setw(8) << setfill(' ') <<"OOC"
				 << left << setw(17)<< setfill(' ') <<"VELY_ERR"
				 << left << setw(8) << setfill(' ') <<"OOC"
				 << left << setw(17)<< setfill(' ') <<"VELZ_ERR"
			     << left << setw(8) << setfill(' ') <<"OOC"
				 << left << setw(17)<< setfill(' ') <<"PRESSURE_ERR"
				 << left << setw(8) << setfill(' ') <<"OOC"<<endl
                 << left << setw(112) << setfill('-') << ""<<endl; 
   
   
   L1_err_file<< left << setw(6) << setfill(' ') << param.Nx_cell[0]
			  << left << setw(6) << setfill(' ') << param.Ny_cell[0]
			  << left << setw(6) << setfill(' ') << param.Nz_cell[0]
			  << scientific
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[0].density
			  << left << setw(8) << setfill(' ') << " - " 
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[0].velocity.x
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[0].velocity.y
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[0].velocity.z
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[0].pressure
			  << left << setw(8) << setfill(' ') << " - "
			  << endl;
   
   
   L2_err_file<< left << setw(6) << setfill(' ') << param.Nx_cell[0]
			  << left << setw(6) << setfill(' ') << param.Ny_cell[0]
			  << left << setw(6) << setfill(' ') << param.Nz_cell[0]
			  << scientific
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[0].density
			  << left << setw(8) << setfill(' ') << " - " 
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[0].velocity.x
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[0].velocity.y
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[0].velocity.z
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[0].pressure
			  << left << setw(8) << setfill(' ') << " - "
			  << endl;
			  
   Linf_err_file<< left << setw(6) << setfill(' ') << param.Nx_cell[0]
			  << left << setw(6) << setfill(' ') << param.Ny_cell[0]
			  << left << setw(6) << setfill(' ') << param.Nz_cell[0]
			  << scientific
			  << left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[0].density
			  << left << setw(8) << setfill(' ') << " - " 
			  << left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[0].velocity.x
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[0].velocity.y
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[0].velocity.z
			  << left << setw(8) << setfill(' ') << " - "
			  << left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[0].pressure
			  << left << setw(8) << setfill(' ') << " - "
			  << endl;			  
   
   double factor,f0=1,f1=1,ratio=1,eps=1.0e-14;
   for(int i=1; i<param.mesh_levels; ++i)
   {
      
      if(param.pdim == 1)
      {
          f0 = pow((param.xmax - param.xmin)/(double)param.Nx_cell[i-1],2.0);
          f1 = pow((param.xmax - param.xmin)/(double)param.Nx_cell[i],2.0); 	    
      }
      else if(param.pdim == 2)
      {
          f0 = pow((param.xmax - param.xmin)/(double)param.Nx_cell[i-1],2.0)
              + pow((param.ymax - param.ymin)/(double)param.Ny_cell[i-1],2.0);
          f1 = pow((param.xmax - param.xmin)/(double)param.Nx_cell[i],2.0)
              + pow((param.ymax - param.ymin)/(double)param.Ny_cell[i],2.0); 	    
      }
      else if(param.pdim == 3)
      {
          f0 = pow((param.xmax - param.xmin)/(double)param.Nx_cell[i-1],2.0)
              + pow((param.ymax - param.ymin)/(double)param.Ny_cell[i-1],2.0)
              + pow((param.zmax - param.zmin)/(double)param.Nz_cell[i-1],2.0);
          f1 = pow((param.xmax - param.xmin)/(double)param.Nx_cell[i],2.0)
              + pow((param.ymax - param.ymin)/(double)param.Ny_cell[i],2.0)
              + pow((param.zmax - param.zmin)/(double)param.Nz_cell[i],2.0); 	    
      }               
      factor = log(sqrt(f1)/sqrt(f0));
      
      // factor = log((double)param.Nx_cell[i-1]*(double)param.Ny_cell[i-1]*(double)param.Nz_cell[i-1]/
//                    ((double)param.Nx_cell[i]*(double)param.Ny_cell[i]*(double)param.Nz_cell[i]));
      
      // L1 Error and OOC
      L1_err_file   << left << setw(6) << setfill(' ') << param.Nx_cell[i]
			        << left << setw(6) << setfill(' ') << param.Ny_cell[i]
			        << left << setw(6) << setfill(' ') << param.Nz_cell[i]
			        << scientific
                    << left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[i].density;
	  if(L1_err[i].density	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L1_err[i].density/L1_err[i-1].density);   			
					
	  L1_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
					<< scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[i].velocity.x;
					
	  if(L1_err[i].velocity.x	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L1_err[i].velocity.x/L1_err[i-1].velocity.x);				
					
	  L1_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[i].velocity.y;
	  
	  if(L1_err[i].velocity.y	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L1_err[i].velocity.y/L1_err[i-1].velocity.y);				
					
	  L1_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[i].velocity.z;
	  
	  if(L1_err[i].velocity.z	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L1_err[i].velocity.z/L1_err[i-1].velocity.z);
	     				
	  L1_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[i].pressure;
	  
	  if(L1_err[i].pressure	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L1_err[i].pressure/L1_err[i-1].pressure);
	     				
	  L1_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
					<< endl;
	  
	  
	  // L2 Error and OOC
      L2_err_file   << left << setw(6) << setfill(' ') << param.Nx_cell[i]
			        << left << setw(6) << setfill(' ') << param.Ny_cell[i]
			        << left << setw(6) << setfill(' ') << param.Nz_cell[i]
			        << scientific
                    << left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[i].density;
	  if(L2_err[i].density	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L2_err[i].density/L2_err[i-1].density);   			
					
	  L2_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[i].velocity.x;
					
	  if(L2_err[i].velocity.x	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L2_err[i].velocity.x/L2_err[i-1].velocity.x);				
					
	  L2_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[i].velocity.y;
	
	  if(L2_err[i].velocity.y	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L2_err[i].velocity.y/L2_err[i-1].velocity.y);
	     
	  L2_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[i].velocity.z;
	  
	  if(L2_err[i].velocity.z	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L2_err[i].velocity.z/L2_err[i-1].velocity.z);   
	     				
	  L2_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[i].pressure;
	  
	  if(L2_err[i].pressure	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L2_err[i].pressure/L2_err[i-1].pressure);
	     				
	  L2_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
					<< endl;
   
	  
	  // Linf Error and OOC
      Linf_err_file   << left << setw(6) << setfill(' ') << param.Nx_cell[i]
			          << left << setw(6) << setfill(' ') << param.Ny_cell[i]
			          << left << setw(6) << setfill(' ') << param.Nz_cell[i]
			          << scientific
                      << left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[i].density;
	  if(Linf_err[i].density	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(Linf_err[i].density/Linf_err[i-1].density);   			
					
	  Linf_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[i].velocity.x;
					
	  if(Linf_err[i].velocity.x	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(Linf_err[i].velocity.x/Linf_err[i-1].velocity.x);				
					
	  Linf_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[i].velocity.y;
	
	  if(Linf_err[i].velocity.y	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(Linf_err[i].velocity.y/Linf_err[i-1].velocity.y);
	     
	  Linf_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[i].velocity.z;
	  
	  if(Linf_err[i].velocity.z	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(Linf_err[i].velocity.z/Linf_err[i-1].velocity.z);   
	     				
	  Linf_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
	                << scientific
					<< left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[i].pressure;
	  
	  if(Linf_err[i].pressure	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(Linf_err[i].pressure/Linf_err[i-1].pressure);
	     				
	  Linf_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
					<< endl;   
   }
                                
}
