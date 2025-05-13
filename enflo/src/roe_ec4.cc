#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// ROE_EC4 flux
//------------------------------------------------------------------------------
void FV::roe_ec4_flux (const vector<PrimVar> &P,
                           const Vector &unit_normal, 
					       Flux &flux) const
{  
   int array_size = P.size();
   assert(array_size>=4);
   Flux flux0;
   vector<PrimVar> P0(2);
   flux.zero();
   
   P0[0] = P[array_size/2-1];
   P0[1] = P[array_size/2];
   roe_ec_flux(P0,unit_normal,flux0);
   flux.sadd(flux0,4.0/3.0);
   
   P0[0] = P[array_size/2-2];
   P0[1] = P[array_size/2];
   roe_ec_flux(P0,unit_normal,flux0);
   flux.sadd(flux0,-1.0/6.0);
   
   P0[0] = P[array_size/2-1];
   P0[1] = P[array_size/2+1];
   roe_ec_flux(P0,unit_normal,flux0);
   flux.sadd(flux0,-1.0/6.0);
}
