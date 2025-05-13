#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// TeCNO scheme for burgers equation
//------------------------------------------------------------------------------
void FV::bur_tecno_flux (const vector<double> &Soln,
                         const Vector &unit_normal,
					     double &flux) const
{  
   int array_size = Soln.size();
   Vector hvec;
   hvec.x = dx, hvec.y = dy, hvec.z = dz;
   double h = hvec*unit_normal;
       
   
   double Ul = Soln[array_size/2-1];
   double Ur = Soln[array_size/2];

   // central flux
   bur_ec_flux(Soln,unit_normal,flux);
                      
   // entropy dissipation
   double D = 0.5*(fabs(Ul) + fabs(Ur));      
							   
   // Reconstructed states
   double Uij,Uji;   
   reconstruct(param.reconstruct_scheme,h,Soln,Uij,Uji); 
   
   flux-=0.5*D*(Uji-Uij);
             
}
