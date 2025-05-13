#include <iostream>
#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;


//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void FV::num_flux(const vector<double> &Soln,
				  const Vector &unit_normal,
				  double &flux) const
{
   double face_area;
   if(unit_normal.x == 1)
     face_area = dy*dz;
   else if(unit_normal.y == 1) 
     face_area = dx*dz;
   else 
     face_area = dx*dy;
   
   switch (param.flux_scheme)
   {
      case lin_ec:
         lin_ec_flux (Soln,unit_normal,flux);
         break;
      case lin_ec4:
         lin_ec4_flux (Soln,unit_normal,flux);
         break;   
      case lin_tecno:
         lin_tecno_flux(Soln,unit_normal,flux);
         break;
      case lin_tecno4:
         lin_tecno4_flux(Soln,unit_normal,flux);
         break;
      case bur_ec:
         bur_ec_flux (Soln,unit_normal,flux);
         break;
      case bur_ec4:
         bur_ec4_flux (Soln,unit_normal,flux);
         break;   
      case bur_tecno:
         bur_tecno_flux(Soln,unit_normal,flux);
         break;
      case bur_tecno4:
         bur_tecno4_flux(Soln,unit_normal,flux);
         break;                    
      default:
         cout << "num_flux: unknown flux " << param.flux_scheme << endl;
         exit (0);
   }
   
   flux*=face_area;
     
}
