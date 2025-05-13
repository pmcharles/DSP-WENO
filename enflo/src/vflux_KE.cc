#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Nodal viscous flux (ES formulation)
// (i,j,k) --> corresponds to node (i+h/2,j+h/2,k+h/2) and will 
// represent the index for the primitive variable at cell-center 
// [ng_cell + i][ng_cell + j][ng_cell +k]
// Algorithm needs one ghost layer around mesh
//------------------------------------------------------------------------------

void FV::nodal_vlux_KE (int i, int j, int k,
                        Flux &xflux, Flux &yflux, Flux &zflux) const
{  
   PrimVar prim_avg;
   Vector Txyz_vel, Tyxz_vel, Tzxy_vel;
   double Txyz_temp, Tyxz_temp, Tzxy_temp;

   assert(i>=-1); assert(i<Nx_cell);
   assert(j>=-1); assert(j<Ny_cell);
   assert(k>=-1); assert(k<Nz_cell);
                
   
   nodal_avg_and_vel_derivatives (i, j, k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
   double T_avg = Temperature(prim_avg); 
   double mu =  param.Viscosity(T_avg);    
   double kappa = mu*param.Cp/param.Pr;
   
   double factx = 1.0/4.0/dx, facty = 1.0/4.0/dy, factz = 1.0/4.0/dz;
   
   // We need derivatives of temperature as well
   
   Txyz_temp =  factx*( Temperature(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1])
                      + Temperature(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k   ]) 
                      + Temperature(primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1])
                      + Temperature(primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k   ]) 
                      - Temperature(primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1])
                      - Temperature(primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k   ]) 
                      - Temperature(primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k +1])
                      - Temperature(primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k   ]) );          
        
   Tyxz_temp =  facty*( Temperature(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1])
                      + Temperature(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k   ]) 
                      + Temperature(primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1])
                      + Temperature(primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k   ]) 
                      - Temperature(primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1])
                      - Temperature(primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k   ]) 
                      - Temperature(primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k +1])
                      - Temperature(primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k   ]) ); 
                      
   Tzxy_temp =  factz*( Temperature(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1])
                      + Temperature(primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1]) 
                      + Temperature(primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1])
                      + Temperature(primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k +1]) 
                      - Temperature(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k  ])
                      - Temperature(primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k  ]) 
                      - Temperature(primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k  ])
                      - Temperature(primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k  ]) );     
   
   double grad_vel = Txyz_vel.x + Tyxz_vel.y + Tzxy_vel.z;       
   double tau_xx = mu*(2*Txyz_vel.x - (2.0/3.0)*grad_vel);
   double tau_yy = mu*(2*Tyxz_vel.y - (2.0/3.0)*grad_vel);
   double tau_zz = mu*(2*Tzxy_vel.z - (2.0/3.0)*grad_vel);
   double tau_xy = mu*(Tyxz_vel.x + Txyz_vel.y);
   double tau_xz = mu*(Txyz_vel.z + Tzxy_vel.x);          
   double tau_yz = mu*(Tyxz_vel.z + Tzxy_vel.y);
   
   //cout << prim_avg.density <<" "<<prim_avg.velocity.x<<" "<<prim_avg.velocity.y<<" "<<prim_avg.velocity.z<<" "<<prim_avg.pressure<<endl;
 

   xflux.mass_flux = 0.0;
   xflux.momentum_flux.x = tau_xx;
   xflux.momentum_flux.y = tau_xy; 
   xflux.momentum_flux.z = tau_xz; 
   xflux.energy_flux  = prim_avg.velocity*xflux.momentum_flux + kappa*Txyz_temp; 

   
   yflux.mass_flux = 0.0;
   yflux.momentum_flux.x = tau_xy;
   yflux.momentum_flux.y = tau_yy; 
   yflux.momentum_flux.z = tau_yz; 
   yflux.energy_flux  = prim_avg.velocity*yflux.momentum_flux + kappa*Tyxz_temp; 

   
   zflux.mass_flux = 0.0;
   zflux.momentum_flux.x = tau_xz;
   zflux.momentum_flux.y = tau_yz; 
   zflux.momentum_flux.z = tau_zz; 
   zflux.energy_flux  = prim_avg.velocity*zflux.momentum_flux + kappa*Tzxy_temp;   
   
   xflux *= -(dy*dz); yflux *= -(dx*dz); zflux *= -(dx*dy); 

}
