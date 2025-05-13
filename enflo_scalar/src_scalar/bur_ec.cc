#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Entropy Conservative Flux for Burgers equation
//------------------------------------------------------------------------------
void FV::bur_ec_flux (const vector<double> &Soln,
                      const Vector &unit_normal,
					  double &flux) const
{  
   int array_size = Soln.size();
   double Ul = Soln[array_size/2-1];
   double Ur = Soln[array_size/2];
   flux = (Ul*Ul + Ur*Ur + Ul*Ur)/6.0;                  
                                                      
}
