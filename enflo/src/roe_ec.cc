#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Roe's Entropy Conservative Flux
//------------------------------------------------------------------------------
void FV::roe_ec_flux (const vector<PrimVar> &P,
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
   
   double bl = sqrt( rhol / pl );
   double br = sqrt( rhor / pr );
   
   // arithmetic avg of parameter vector
   double z1b = 0.5*(bl + br);
   Vector zvb;
   zvb.equ(Pl.velocity,Pr.velocity,0.5*bl,0.5*br);
   double z5b = 0.5*(bl*pl + br*pr);
   double zvb_normal = zvb * unit_normal;
   
   //logarithmic avg of parameter vector
   double z1l = logavg( bl, br );
   double z5l = logavg( bl*pl, br*pr );
   
   // double rho = z1b*z5l;
//    Vector vel;
//    vel.equ(zvb,1.0/z1b);
//    double vel_normal = vel*unit_normal;
//    double p1 = z5b/z1b;
//    double p2 = ((param.GAMMA+1.0)*z5l/z1l + (param.GAMMA-1)*z5b/z1b)/2.0/param.GAMMA;
//    double a = sqrt(param.GAMMA*p2/rho);
//    double H = a*a/(param.GAMMA-1) + vel.square()*0.5;
   
   
   // central flux
   // flux.mass_flux = rho*vel_normal;
//    flux.momentum_flux.equ(unit_normal,vel,p1,flux.mass_flux);
//    flux.energy_flux = rho*H*vel_normal;
   
   // central flux
   flux.mass_flux = z5l * zvb_normal;
   flux.momentum_flux.equ(unit_normal,zvb,z5b/z1b,flux.mass_flux/z1b);
   flux.energy_flux = 0.5*( (param.GAMMA+1.0)/(param.GAMMA-1.0)* flux.mass_flux/z1l + 
                      flux.momentum_flux * zvb)/z1b;  
                      
   flux*=face_area;                   
                                                      
}
