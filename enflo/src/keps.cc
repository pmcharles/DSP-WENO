#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// KEPS flux (kinetic energy and approximately entropy preserving flux)
//------------------------------------------------------------------------------
void FV::keps_flux (const vector<PrimVar> &P,
                           const Vector &unit_normal, 
					       Flux &flux) const
{  
   int s = P.size();
   
   
   PrimVar Pl = P[s/2-1];
   PrimVar Pr = P[s/2];
   double rhol = Pl.density;
   double rhor = Pr.density;
   double rho = 0.5*(rhol+rhor);
   double templ = Pl.pressure/(param.gas_const * rhol);
   double tempr = Pr.pressure/(param.gas_const * rhor);
   double vel2= 0.5 * (Pl.velocity.square() + Pr.velocity.square());
   double betal = 0.5 / (param.gas_const *templ);
   double betar = 0.5 / (param.gas_const *tempr);
   double beta = 0.5*(betal+betar);
   double p     = 0.5 * rho / beta;
   Vector vel;
   vel.equ(Pl.velocity,Pr.velocity,0.5,0.5);
   double vel_normal = vel * unit_normal;

   // central flux
   flux.mass_flux = rho * vel_normal;
   flux.momentum_flux.equ(unit_normal,vel,p,flux.mass_flux);
   flux.energy_flux = 0.5 * ( 1.0/((param.GAMMA-1.0)*beta) - vel2) * flux.mass_flux + 
                      flux.momentum_flux *  vel;   
    
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
