#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Entropy Conservative Flux for Linear Advection
//------------------------------------------------------------------------------
void FV::lin_ec_flux (const vector<double> &Soln,
                      const Vector &unit_normal,
					  double &flux) const
{  
   int array_size = Soln.size();
   double Ul = Soln[array_size/2-1];
   double Ur = Soln[array_size/2];
   double vel = unit_normal*param.lin_vel;
   flux = 0.5*vel*(Ul + Ur);                  
                                                      
}
