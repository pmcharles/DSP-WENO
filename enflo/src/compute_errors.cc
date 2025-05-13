#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Computing errors
//------------------------------------------------------------------------------
void FV::compute_error(const int &ref_level)
{
   PrimVar abs_diff,exact_pw;
   L1_err[ref_level]   = 0.0;
   L2_err[ref_level]   = 0.0;
   Linf_err[ref_level] = 0.0;
   
   for(int i=0; i<Nx_cell; ++i)
      for(int j=0; j<Ny_cell; ++j)
		 for(int k=0; k<Nz_cell; ++k)
		 {
			param.exact_soln.value (cell_cc[i][j][k],time,exact_pw);
			abs_diff.adiff_eq(exact_pw,primitive[i+ng_cell][j+ng_cell][k+ng_cell]);
			L1_err[ref_level] += abs_diff;
			L2_err[ref_level].sqadd(abs_diff);
			Linf_err[ref_level].max(abs_diff);
		 }
   L1_err[ref_level] *= dx*dy*dz;
   L2_err[ref_level] *= dx*dy*dz;
      
   L2_err[ref_level].density    = sqrt(L2_err[ref_level].density);
   L2_err[ref_level].velocity.x = sqrt(L2_err[ref_level].velocity.x);
   L2_err[ref_level].velocity.y = sqrt(L2_err[ref_level].velocity.y); 
   L2_err[ref_level].velocity.z = sqrt(L2_err[ref_level].velocity.z); 
   L2_err[ref_level].pressure = sqrt(L2_err[ref_level].pressure);                          
}