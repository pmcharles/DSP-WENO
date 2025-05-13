#include <iostream>
#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;


//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void FV::num_flux(const vector<PrimVar> &P,
				  const vector<EntVar> &V,
				  const Vector &unit_normal,
				  Flux &flux) const
{
   switch (param.flux_scheme)
   {
      case kepes_tecno_roe:
         kepes_tecno_roe_flux (P,V,unit_normal,flux);
         break;
      case kepes_tecno_rusanov:
         kepes_tecno_rusanov_flux (P,V,unit_normal,flux);
         break;   
      case kepec:
         kepec_flux(P,unit_normal,flux);
         break;
      case kepec4:
         kepec4_flux(P,unit_normal,flux);
         break;
      case kepes_tecno4_roe:
         kepes_tecno4_roe_flux(P,V,unit_normal,flux);
         break;
      case kepes_tecno4_rusanov:
         kepes_tecno4_rusanov_flux(P,V,unit_normal,flux);
         break;
      case kep:
         kep_flux(P,unit_normal,flux);
         break;
      case kep4:
         kep4_flux(P,unit_normal,flux);
         break; 
      case keps:
         keps_flux(P,unit_normal,flux);
         break;
      case keps4:
         keps4_flux(P,unit_normal,flux);
         break;     
      case kep_tecno_roe:
         kep_tecno_roe_flux (P,V,unit_normal,flux);
         break;  
      case roe:
         roe_flux(P,unit_normal,flux);
         break;     
      case roe_ec:
         roe_ec_flux(P,unit_normal,flux);
         break;               
      case roe_ec4:
         roe_ec4_flux(P,unit_normal,flux);
         break;  
      case roe_tecno_roe:
         roe_tecno_roe_flux(P,V,unit_normal,flux);
         break;               
      case roe_tecno4_roe:
         roe_tecno4_roe_flux(P,V,unit_normal,flux);
         break;                     
      default:
         cout << "num_flux: unknown flux " << param.flux_scheme << endl;
         exit (0);
   }
     
}
