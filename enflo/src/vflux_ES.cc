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

void FV::nodal_vlux_ES (int i, int j, int k,
                 Flux &xflux, Flux &yflux, Flux &zflux) const
{  
   EntVar ent_avg, Txyz_ent, Tyxz_ent, Tzxy_ent;

   assert(i>=-1); assert(i<Nx_cell);
   assert(j>=-1); assert(j<Ny_cell);
   assert(k>=-1); assert(k<Nz_cell);
           
   double fact1 = 1.0/8.0;        

   ent_avg.equ(entropy_var[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1],
               entropy_var[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1], 
               entropy_var[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1],
               entropy_var[ng_cell + i +1][ng_cell + j +1][ng_cell +k   ],
               fact1, fact1, fact1, fact1);
   ent_avg.sadd(entropy_var[ng_cell + i   ][ng_cell + j   ][ng_cell +k +1],
                entropy_var[ng_cell + i   ][ng_cell + j +1][ng_cell +k   ],
                entropy_var[ng_cell + i +1][ng_cell + j   ][ng_cell +k   ],
                entropy_var[ng_cell + i   ][ng_cell + j   ][ng_cell +k   ],
                fact1, fact1, fact1, fact1);    
   
   double factx = 1.0/4.0/dx, facty = 1.0/4.0/dy, factz = 1.0/4.0/dz;
   

   Txyz_ent.equ(entropy_var[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1] ,
                entropy_var[ng_cell + i +1][ng_cell + j +1][ng_cell +k   ] , 
                entropy_var[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1] ,
                entropy_var[ng_cell + i +1][ng_cell + j   ][ng_cell +k   ] ,
                factx, factx, factx, factx); 
   Txyz_ent.sadd(entropy_var[ng_cell + i][ng_cell + j +1][ng_cell +k +1] ,
                 entropy_var[ng_cell + i][ng_cell + j +1][ng_cell +k   ] , 
                 entropy_var[ng_cell + i][ng_cell + j   ][ng_cell +k +1] ,
                 entropy_var[ng_cell + i][ng_cell + j   ][ng_cell +k   ] ,
                 -factx, -factx, -factx, -factx);
        
              
   Tyxz_ent.equ(entropy_var[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1] ,
                entropy_var[ng_cell + i +1][ng_cell + j +1][ng_cell +k   ] , 
                entropy_var[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1] ,
                entropy_var[ng_cell + i   ][ng_cell + j +1][ng_cell +k   ] ,
                facty, facty, facty, facty);
   Tyxz_ent.sadd(entropy_var[ng_cell + i +1][ng_cell + j][ng_cell +k +1] ,
                 entropy_var[ng_cell + i +1][ng_cell + j][ng_cell +k   ] , 
                 entropy_var[ng_cell + i   ][ng_cell + j][ng_cell +k +1] ,
                 entropy_var[ng_cell + i   ][ng_cell + j][ng_cell +k   ] ,
                 -facty, -facty, -facty, -facty);
                      
   Tzxy_ent.equ(entropy_var[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1] ,
                entropy_var[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1] , 
                entropy_var[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1] ,
                entropy_var[ng_cell + i   ][ng_cell + j   ][ng_cell +k +1] ,
                factz, factz, factz, factz);
   Tzxy_ent.sadd(entropy_var[ng_cell + i +1][ng_cell + j +1][ng_cell +k] ,
                 entropy_var[ng_cell + i +1][ng_cell + j   ][ng_cell +k] , 
                 entropy_var[ng_cell + i   ][ng_cell + j +1][ng_cell +k] ,
                 entropy_var[ng_cell + i   ][ng_cell + j   ][ng_cell +k] ,
                 -factz, -factz, -factz, -factz);
   
   double T_avg = -1/ent_avg.entl/param.gas_const;
   double mu =  param.Viscosity(T_avg);    
   double kappa = mu*param.Cp/param.Pr;
   double c1 = 4.0*mu/3.0, c2 = -2.0*mu/3.0, c3 = kappa/param.gas_const, e4 = 1.0/ent_avg.entl, 
          e1 = ent_avg.entv.x*e4, e2 = ent_avg.entv.y*e4, e3 = ent_avg.entv.z*e4;
   
   // X-FLUX 
   // K11*V_x
   xflux.mass_flux = 0.0;
   xflux.momentum_flux.x = c1*e4*(-Txyz_ent.entv.x + e1*Txyz_ent.entl);
   xflux.momentum_flux.y = mu*e4*(-Txyz_ent.entv.y + e2*Txyz_ent.entl);
   xflux.momentum_flux.z = mu*e4*(-Txyz_ent.entv.z + e3*Txyz_ent.entl);
   xflux.energy_flux     = e4*(c1*e1*Txyz_ent.entv.x + mu*e2*Txyz_ent.entv.y + mu*e3*Txyz_ent.entv.z  
                            +(-c1*e1*e1 - mu*(e2*e2 + e3*e3) + c3*e4)*Txyz_ent.entl);
                            
   // K12*V_y
   xflux.momentum_flux.x += c2*e4*(-Tyxz_ent.entv.y + e2*Tyxz_ent.entl);
   xflux.momentum_flux.y += mu*e4*(-Tyxz_ent.entv.x + e1*Tyxz_ent.entl);
   xflux.momentum_flux.z += 0;
   xflux.energy_flux     += e4*(mu*e2*Tyxz_ent.entv.x + c2*e1*Tyxz_ent.entv.y - (mu+c2)*e1*e2*Tyxz_ent.entl);
   
   // K13*V_z
   xflux.momentum_flux.x += c2*e4*(-Tzxy_ent.entv.z + e3*Tzxy_ent.entl);
   xflux.momentum_flux.y += 0;
   xflux.momentum_flux.z += mu*e4*(-Tzxy_ent.entv.x + e1*Tzxy_ent.entl);
   xflux.energy_flux     += e4*(mu*e3*Tzxy_ent.entv.x + c2*e1*Tzxy_ent.entv.z - (mu+c2)*e1*e3*Tzxy_ent.entl);    
   
   
   
   // Y-FLUX 
   // K21*V_x
   yflux.momentum_flux.x += mu*e4*(-Txyz_ent.entv.y + e2*Txyz_ent.entl);
   yflux.momentum_flux.y += c2*e4*(-Txyz_ent.entv.x + e1*Txyz_ent.entl);
   yflux.momentum_flux.z += 0;
   yflux.energy_flux     += e4*(c2*e2*Txyz_ent.entv.x + mu*e1*Txyz_ent.entv.y - (mu+c2)*e1*e2*Txyz_ent.entl);
   
                            
   // K22*V_y
   yflux.mass_flux = 0.0;
   yflux.momentum_flux.x = mu*e4*(-Tyxz_ent.entv.x + e1*Tyxz_ent.entl);
   yflux.momentum_flux.y = c1*e4*(-Tyxz_ent.entv.y + e2*Tyxz_ent.entl);
   yflux.momentum_flux.z = mu*e4*(-Tyxz_ent.entv.z + e3*Tyxz_ent.entl);
   yflux.energy_flux     = e4*(mu*e1*Tyxz_ent.entv.x + c1*e2*Tyxz_ent.entv.y + mu*e3*Tyxz_ent.entv.z  
                            +(-c1*e2*e2 - mu*(e1*e1 + e3*e3) + c3*e4)*Tyxz_ent.entl);
   
   // K23*V_z
   yflux.momentum_flux.x += 0;
   yflux.momentum_flux.y += c2*e4*(-Tzxy_ent.entv.z + e3*Tzxy_ent.entl);
   yflux.momentum_flux.z += mu*e4*(-Tzxy_ent.entv.y + e2*Tzxy_ent.entl);
   yflux.energy_flux     += e4*(mu*e3*Tzxy_ent.entv.y + c2*e2*Tzxy_ent.entv.z - (mu+c2)*e2*e3*Tzxy_ent.entl); 
   
   
   // Z-FLUX 
   // K31*V_x
   zflux.momentum_flux.x += mu*e4*(-Txyz_ent.entv.z + e3*Txyz_ent.entl);
   zflux.momentum_flux.y += 0;
   zflux.momentum_flux.z += c2*e4*(-Txyz_ent.entv.x + e1*Txyz_ent.entl);
   zflux.energy_flux     += e4*(c2*e3*Txyz_ent.entv.x + mu*e1*Txyz_ent.entv.z - (mu+c2)*e1*e3*Txyz_ent.entl);
   
                            
   // K32*V_y
   zflux.momentum_flux.x += 0;
   zflux.momentum_flux.y += mu*e4*(-Tyxz_ent.entv.z + e3*Tyxz_ent.entl);
   zflux.momentum_flux.z += c2*e4*(-Tyxz_ent.entv.y + e2*Tyxz_ent.entl);
   zflux.energy_flux     += e4*(c2*e3*Tyxz_ent.entv.y + mu*e2*Tyxz_ent.entv.z - (mu+c2)*e2*e3*Tyxz_ent.entl); 
   
   
   // K33*V_z
   zflux.mass_flux = 0.0;
   zflux.momentum_flux.x = mu*e4*(-Tzxy_ent.entv.x + e1*Tzxy_ent.entl);
   zflux.momentum_flux.y = mu*e4*(-Tzxy_ent.entv.y + e2*Tzxy_ent.entl);
   zflux.momentum_flux.z = c1*e4*(-Tzxy_ent.entv.z + e3*Tzxy_ent.entl);
   zflux.energy_flux     = e4*(mu*e1*Tzxy_ent.entv.x + mu*e2*Tzxy_ent.entv.y + c1*e3*Tzxy_ent.entv.z  
                            +(-c1*e3*e3 - mu*(e1*e1 + e2*e2) + c3*e4)*Tzxy_ent.entl);
   

   xflux *= -(dy*dz); yflux *= -(dx*dz); zflux *= -(dx*dy); 
}
