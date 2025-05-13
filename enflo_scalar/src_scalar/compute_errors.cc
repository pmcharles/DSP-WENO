#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Computing errors
//------------------------------------------------------------------------------
void FV::compute_error(const int &ref_level)
{
   double abs_diff,exact_pw;
   L1_err[ref_level]   = 0.0;
   L2_err[ref_level]   = 0.0;
   Linf_err[ref_level] = 0.0;
   
   for(int i=0; i<Nx_cell; ++i)
      for(int j=0; j<Ny_cell; ++j)
		 for(int k=0; k<Nz_cell; ++k)
		 {
			param.exact_soln.value (cell_cc[i][j][k],time,exact_pw);
			abs_diff = fabs(exact_pw - solution[i+ng_cell][j+ng_cell][k+ng_cell]);
			L1_err[ref_level] += abs_diff;
			L2_err[ref_level] += pow(abs_diff,2);
			Linf_err[ref_level] = max(Linf_err[ref_level],abs_diff);
		 }
   L1_err[ref_level] *= dx*dy*dz;
   L2_err[ref_level] *= dx*dy*dz;
      
   L2_err[ref_level] = sqrt(L2_err[ref_level]);                          
}