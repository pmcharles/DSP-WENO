#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Kinetic Energy preserving flux of Jameson flux
//------------------------------------------------------------------------------
void FV::kep_flux (const vector<PrimVar> &P,
                           const Vector &unit_normal, 
					       Flux &flux) const
{  
   int s = P.size();
   
   PrimVar Pl = P[s/2-1];
   PrimVar Pr = P[s/2];
   
   double rhol = Pl.density;
   double rhor = Pr.density;
   double rho  = 0.5*(rhol+rhor);
   double p     = 0.5 * (Pl.pressure + Pr.pressure);
   Vector vel;
   vel.equ(Pl.velocity,Pr.velocity,0.5,0.5);
   double H = 0.5*(Enthalpy(Pl) + Enthalpy(Pr));
   
   double vel_normal = vel * unit_normal;
   

   // central flux
   flux.mass_flux = rho * vel_normal;
   flux.momentum_flux.equ(unit_normal,vel,p,flux.mass_flux);
   flux.energy_flux = rho*H*vel_normal;   
    
   double face_area;
   if(unit_normal.x == 1)
     face_area = dy*dz;
   else if(unit_normal.y == 1) 
     face_area = dx*dz;
   else 
     face_area = dx*dy;    
   flux*=face_area;                   
    //cout<<flux.mass_flux<<" "<<flux.momentum_flux.x<<" "<<flux.momentum_flux.y<<" "<<flux.energy_flux<<endl;                                   
}
