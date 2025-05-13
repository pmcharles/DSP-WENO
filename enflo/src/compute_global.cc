#include <cmath>
#include <iomanip>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Computing and saving global quantities
//------------------------------------------------------------------------------
void FV::compute_global()
{

   
   // Cell-center quantities
   // KE
   // total_ent
   // KE_diss_pr (already evaluated)
   KE = 0.0;
   total_ent = 0.0;
   for(int i=0; i<Nx_cell; ++i)
      for(int j=0; j<Ny_cell; ++j)
         for(int k=0; k<Nz_cell; ++k)
		 {
			double s = Entropy (primitive[i+ng_cell][j+ng_cell][k+ng_cell]);
			total_ent += -primitive[i+ng_cell][j+ng_cell][k+ng_cell].density*s/(param.GAMMA -1.0); 
			KE += primitive[i+ng_cell][j+ng_cell][k+ng_cell].density*primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity.square();

		 }
         
   // Evaluated at nodes
   // KE_diss_stress
   // Enstrophy
   Enstrophy = 0.0;
      
   // Loop over nodes
   double fact1 = 1.0/8.0;   
   double den_avg;
   for(int i=0; i<Nx_cell+1; ++i)
      for(int j=0; j<Ny_cell+1; ++j)
         for(int k=0; k<Nz_cell+1; ++k)
		 {
			
            den_avg = fact1*(  primitive[ng_cell + i -1][ng_cell + j -1][ng_cell +k -1].density
                             + primitive[ng_cell + i   ][ng_cell + j -1][ng_cell +k -1].density 
                             + primitive[ng_cell + i -1][ng_cell + j   ][ng_cell +k -1].density
                             + primitive[ng_cell + i -1][ng_cell + j -1][ng_cell +k   ].density
                             + primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k -1].density
                             + primitive[ng_cell + i   ][ng_cell + j -1][ng_cell +k   ].density
                             + primitive[ng_cell + i -1][ng_cell + j   ][ng_cell +k   ].density
                             + primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k   ].density
                             );

            Enstrophy += den_avg*vorticity[i][j][k].square()*node_w[i][j][k];
            
            
         }   
         
   double cell_volume = dx*dy*dz;
   double nodal_volume = (param.xmax - param.xmin) / (Nx_cell+1) *
                         (param.ymax - param.ymin) / (Ny_cell+1) *
                         (param.zmax - param.zmin) / (Nz_cell+1);
         
   total_ent *= cell_volume;
   KE *= 0.5*cell_volume;
   Enstrophy *=0.5*nodal_volume;
   // KE_diss rate have already been mulitplied with the volume
   
   
   if(time == 0)
   {
      total_ent0 = total_ent;
	  if(total_ent0 == 0)
	  {
		 cout<<"WARNING: Initial total entropy is zero!! Re-setting it to 1.0"<<endl;
		 total_ent0 = 1.0;
	  }  
   }
   
      
   glob_file   << time << " "
               << setprecision(10)
               << total_ent << " "
               << (total_ent - total_ent0)/fabs(total_ent0)<< " "
               << KE << " "
               << KE_diss_pr << " "
               << KE_diss_stress << " "
               << Enstrophy << " "
               << endl;  
}
