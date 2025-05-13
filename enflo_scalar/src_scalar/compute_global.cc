#include <cmath>
#include <iomanip>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Computing and saving global quantities
//------------------------------------------------------------------------------
void FV::compute_global()
{   
   total_ent = 0.0;
   for(int i=0; i<Nx_cell; ++i)
      for(int j=0; j<Ny_cell; ++j)
         for(int k=0; k<Nz_cell; ++k)
		 {
			total_ent += 0.5*pow(solution[i+ng_cell][j+ng_cell][k+ng_cell],2); 
		 }
         
   double cell_volume = dx*dy*dz;
   total_ent *= cell_volume;   
   
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
               << (total_ent - total_ent0)/total_ent0<< " "
               << endl;  
}
