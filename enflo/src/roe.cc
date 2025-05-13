#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Roe's Original Flux (witout entropy fix)
//------------------------------------------------------------------------------
void FV::roe_flux (const vector<PrimVar> &P,
                           const Vector &unit_normal,
					       Flux &flux) const
{  
   int array_size = P.size();
   double face_area;
   if(unit_normal.x == 1)
     face_area = dy*dz;
   else if(unit_normal.y == 1) 
     face_area = dx*dz;
   else 
     face_area = dx*dy; 
   
   
   
   PrimVar Pl = P[array_size/2-1];
   PrimVar Pr = P[array_size/2];
   double rhol = Pl.density;
   double rhor = Pr.density;
   double pl = Pl.pressure;
   double pr = Pr.pressure;
   double hl = param.GAMMA*pl/rhol/(param.GAMMA-1.0) + 0.5*Pl.velocity.square();
   double hr = param.GAMMA*pr/rhor/(param.GAMMA-1.0) + 0.5*Pr.velocity.square();
   
   // Roe average
   double rhol_sqrt = sqrt(rhol);
   double rhor_sqrt = sqrt(rhor);
   double fact_l    = rhol_sqrt/(rhol_sqrt + rhor_sqrt);
   double fact_r    = 1.0 - fact_l;
   double rho       = rhol_sqrt*rhor_sqrt;
   Vector vel;
   vel.equ(Pl.velocity,Pr.velocity,fact_l,fact_r);
   double h = hl*fact_l + hr*fact_r;
   double vel_normal = vel*unit_normal;
   double a = sqrt( (param.GAMMA - 1.0)*(h - 0.5*vel.square()) );
   
   double dp = pr -pl;
   double vel_l_normal = Pl.velocity*unit_normal;
   double vel_r_normal = Pr.velocity*unit_normal;
   double dV = vel_r_normal - vel_l_normal;
   
   if(vel_normal >=0.0)
   {
       double lambda = vel_normal - a;
       double coeff  = 0.5*(dp - rho*a*dV)/(a*a);
       double factor = min(lambda,0.0)*coeff;
       
       //Left flux
       flux.mass_flux = rhol*vel_l_normal;
       flux.momentum_flux.equ(unit_normal,Pl.velocity,pl,rhol*vel_l_normal);
       flux.energy_flux = hl*flux.mass_flux;
       
       // Upwinding 
       flux.mass_flux += factor;
       flux.momentum_flux.sadd(vel,unit_normal,factor,-a*factor);
       flux.energy_flux += (h - a*vel_normal)*factor;
   }
   else
   {
       double lambda = vel_normal + a;
       double coeff  = 0.5*(dp + rho*a*dV)/(a*a);
       double factor = max(lambda,0.0)*coeff;
       
       //Right flux
       flux.mass_flux = rhor*vel_r_normal;
       flux.momentum_flux.equ(unit_normal,Pr.velocity,pr,rhor*vel_r_normal);
       flux.energy_flux = hr*flux.mass_flux;
       
       // Upwinding 
       flux.mass_flux -= factor;
       flux.momentum_flux.sadd(vel,unit_normal,-factor,-a*factor);
       flux.energy_flux -= (h + a*vel_normal)*factor;
   }
                      
   flux*=face_area;                   
                                                      
}
