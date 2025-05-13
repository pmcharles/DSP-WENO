#include <cmath>
#include <iomanip>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Computing OOC
//------------------------------------------------------------------------------
void FV::compute_OOC()
{
   
   L1_err_file << left << setw(40) << setfill('-') <<""<<endl;
   L1_err_file << left << setw(6) << setfill(' ') << "Nx"
               << left << setw(6) << setfill(' ') << "Ny" 
               << left << setw(6) << setfill(' ') << "Nz" 
			   << left << setw(17)<< setfill(' ') <<"ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"<<endl
               << left << setw(40) << setfill('-') << ""<<endl; 
   
   L2_err_file << left << setw(40) << setfill('-') <<""<<endl;
   L2_err_file << left << setw(6) << setfill(' ') << "Nx"
               << left << setw(6) << setfill(' ') << "Ny" 
               << left << setw(6) << setfill(' ') << "Nz"
			   << left << setw(17)<< setfill(' ') <<"ERR"
			   << left << setw(8) << setfill(' ') <<"OOC"<<endl
               << left << setw(40) << setfill('-') << ""<<endl; 
   
   Linf_err_file << left << setw(40) << setfill('-') <<""<<endl;
   Linf_err_file << left << setw(6) << setfill(' ') << "Nx"
				 << left << setw(6) << setfill(' ') << "Ny" 
				 << left << setw(6) << setfill(' ') << "Nz"
				 << left << setw(17)<< setfill(' ') <<"ERR"
				 << left << setw(8) << setfill(' ') <<"OOC"<<endl
                 << left << setw(40) << setfill('-') << ""<<endl; 
   
   
   L1_err_file<< left << setw(6) << setfill(' ') << param.Nx_cell[0]
			  << left << setw(6) << setfill(' ') << param.Ny_cell[0]
			  << left << setw(6) << setfill(' ') << param.Nz_cell[0]
			  << scientific
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[0]
			  << left << setw(8) << setfill(' ') << " - "
			  << endl;
   
   
   L2_err_file<< left << setw(6) << setfill(' ') << param.Nx_cell[0]
			  << left << setw(6) << setfill(' ') << param.Ny_cell[0]
			  << left << setw(6) << setfill(' ') << param.Nz_cell[0]
			  << scientific
			  << left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[0]
			  << left << setw(8) << setfill(' ') << " - "
			  << endl;
			  
   Linf_err_file<< left << setw(6) << setfill(' ') << param.Nx_cell[0]
			  << left << setw(6) << setfill(' ') << param.Ny_cell[0]
			  << left << setw(6) << setfill(' ') << param.Nz_cell[0]
			  << scientific
			  << left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[0]
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
                    << left << setw(17) << setfill(' ') << setprecision(5) <<L1_err[i];
	  if(L1_err[i]	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L1_err[i]/L1_err[i-1]);   			
					
	  L1_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
					<< endl;
	  
	  
	  // L2 Error and OOC
      L2_err_file   << left << setw(6) << setfill(' ') << param.Nx_cell[i]
			        << left << setw(6) << setfill(' ') << param.Ny_cell[i]
			        << left << setw(6) << setfill(' ') << param.Nz_cell[i]
			        << scientific
                    << left << setw(17) << setfill(' ') << setprecision(5) <<L2_err[i];
	  if(L2_err[i]	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(L2_err[i]/L2_err[i-1]);   			
					
	  L2_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
					<< endl;
   
	  
	  // Linf Error and OOC
      Linf_err_file   << left << setw(6) << setfill(' ') << param.Nx_cell[i]
			          << left << setw(6) << setfill(' ') << param.Ny_cell[i]
			          << left << setw(6) << setfill(' ') << param.Nz_cell[i]
			          << scientific
                      << left << setw(17) << setfill(' ') << setprecision(5) <<Linf_err[i];
	  if(Linf_err[i]	< eps)
	     ratio = 0.0;
	  else
	     ratio = log(Linf_err[i]/Linf_err[i-1]);   			
					
	  Linf_err_file	<< fixed
	                << left << setw(8) << setfill(' ') << setprecision(4) <<ratio/factor
					<< endl;   
   }
                                
}
