#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// EC4 flux for linear advection
//------------------------------------------------------------------------------
void FV::lin_ec4_flux (const vector<double> &Soln,
                       const Vector &unit_normal, 
					   double &flux) const
{  
   int array_size = Soln.size();
   assert(array_size>=4);
   double flux0;
   vector<double> U0(2);
   flux = 0.0;
   
   U0[0] = Soln[array_size/2-1];
   U0[1] = Soln[array_size/2];
   lin_ec_flux(U0,unit_normal,flux0);
   flux+=flux0*4.0/3.0;
   
   U0[0] = Soln[array_size/2-2];
   U0[1] = Soln[array_size/2];
   lin_ec_flux(U0,unit_normal,flux0);
   flux+=flux0*(-1.0/6.0);
   
   U0[0] = Soln[array_size/2-1];
   U0[1] = Soln[array_size/2+1];
   lin_ec_flux(U0,unit_normal,flux0);
   flux+=flux0*(-1.0/6.0);
}
