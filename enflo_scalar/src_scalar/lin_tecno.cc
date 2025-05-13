#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// TeCNO scheme for linear advection
//------------------------------------------------------------------------------
void FV::lin_tecno_flux (const vector<double> &Soln,
                         const Vector &unit_normal,
					     double &flux) const
{  
   Vector hvec;
   hvec.x = dx, hvec.y = dy, hvec.z = dz;
   double h = hvec*unit_normal;

   // central flux
   lin_ec_flux(Soln,unit_normal,flux);
                      
   // entropy dissipation
   double D = fabs(unit_normal*param.lin_vel);      
							   
   // Reconstructed states
   double Uij,Uji;   
   reconstruct(param.reconstruct_scheme,h,Soln,Uij,Uji); 
   
   flux-=0.5*D*(Uji-Uij);
             
}
