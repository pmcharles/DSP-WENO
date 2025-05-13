#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// ROE_EC TeCNO4 ROE Flux
//------------------------------------------------------------------------------
void FV::roe_tecno4_roe_flux (const vector<PrimVar> &P,
                           const vector<EntVar> &V,
                           const Vector &unit_normal,
					       Flux &flux) const
{  
   int array_size = P.size();
   double face_area;
   double h;
   if(unit_normal.x == 1)
   {
     face_area = dy*dz;
     h = dx;
   }  
   else if(unit_normal.y == 1) 
   {
     face_area = dx*dz;
     h = dy;
   } 
   else 
   {
     face_area = dx*dy; 
     h = dz;
   } 
       
   
   PrimVar Pl = P[array_size/2-1];
   PrimVar Pr = P[array_size/2];
   double rhol = Pl.density;
   double rhor = Pr.density;
   double rho = logavg( rhol, rhor );
   double templ = Pl.pressure/(param.gas_const * rhol);
   double tempr = Pr.pressure/(param.gas_const * rhor);
   //double vel2= 0.5 * (Pl.velocity.square() + Pr.velocity.square());
   double betal = 0.5 / (param.gas_const *templ);
   double betar = 0.5 / (param.gas_const *tempr);
   double beta = logavg(betal, betar);
   double p     = 0.5 * (rhol + rhor) / (betal + betar);
   Vector vel;
   vel.equ(Pl.velocity,Pr.velocity,0.5,0.5);
   double vel_normal = vel * unit_normal;

   // central flux
   roe_ec4_flux(P,unit_normal,flux);
   flux*=(1/face_area);
                      
   // entropy dissipation
   // eigenvectors
   double a   = sqrt(0.5 * param.GAMMA / beta);
   double H = a*a/(param.GAMMA-1.0) + 0.5*vel.square();
   double v1 = vel.x * unit_normal.y - vel.y * unit_normal.x;
   double v2 = vel.z * unit_normal.x - vel.x * unit_normal.z;

   double R[][NVAR] = {
		  {                1         ,          1         ,           0         ,          0       ,             1           },
		  {vel.x - a*unit_normal.x   ,         vel.x      ,    unit_normal.y    , -unit_normal.z   , vel.x + a*unit_normal.x }, 
		  {vel.y - a*unit_normal.y   ,         vel.y      ,    -unit_normal.x   ,          0       , vel.y + a*unit_normal.y },
		  {vel.z - a*unit_normal.z   ,         vel.z      ,           0         ,    unit_normal.x , vel.z + a*unit_normal.z },
		  {H     - a*vel_normal      ,    0.5*vel.square(),          v1         ,          v2      ,  H     + a*vel_normal   } 
	   };

   double const1 = sqrt(rho/param.GAMMA);
   double psqrt = sqrt(p);
   double S_sqrt[] = { sqrt(0.5)*const1, sqrt(param.GAMMA-1.0)*const1, psqrt, psqrt, sqrt(0.5)*const1 };
   for(int i=0; i<NVAR; ++i)
	  for(int j=0; j<NVAR; ++j)
		 R[j][i] *= S_sqrt[i];

   // eigenvalues
   double vnl = Pl.velocity  * unit_normal;
   double vnr = Pr.velocity * unit_normal;
   double al  = sqrt(0.5 * param.GAMMA / betal);
   double ar  = sqrt(0.5 * param.GAMMA / betar);

   // ROE
   double LambdaL[] = { vnl - al, vnl, vnl, vnl, vnl + al };
   double LambdaR[] = { vnr - ar, vnr, vnr, vnr, vnr + ar };
   double l2, l3;
   l2 = l3 = fabs(vel_normal);
   
   double ev_mod_beta =0.0;
   double ev_mod_alpha = 0.0;
   
   if(param.EC1_fix)
      ev_mod_alpha = 1.0/6.0;

   double Lambda[]  = { (1+ev_mod_beta)*fabs(vel_normal - a) + ev_mod_alpha*fabs(LambdaL[0]-LambdaR[0]), 
						l2,
						l3,
						l3,
						(1+ev_mod_beta)*fabs(vel_normal + a) + ev_mod_alpha*fabs(LambdaL[3]-LambdaR[3])};

   double Diff[NVAR];             
							   
   // Scaled entropy variables
   vector< vector<double> > Z;
   Z.resize(array_size);
   for(int i=0; i<array_size; ++i)
   {	 
      Z[i].resize(NVAR,0.0);
   }
   
   for(int k=0; k<array_size; ++k)
   {
	  for(int i=0; i<NVAR; ++i)
	  {
		   Z[k][i] += R[0][i]*V[k].entf 
		             +R[1][i]*V[k].entv.x
		             +R[2][i]*V[k].entv.y
		             +R[3][i]*V[k].entv.z
		             +R[4][i]*V[k].entl;
	  }     
   }	
   

   // Reconstructed states
   double Zij[NVAR],Zji[NVAR];   
   reconstruct(param.reconstruct_scheme,h,Z,Zij,Zji); 

   // Jump in reconstructed states
   double dZij[NVAR];
   for(int i=0; i<NVAR; ++i)
   {   
	  dZij[i] = Zji[i] - Zij[i] ; 
	  Diff[i] = 0;
   }    
   
   
   

   // diffusive flux = R * Lambda * dZij
   for(int i=0; i<NVAR; ++i)
   {
	 for(int j=0; j<NVAR; ++j)
		 Diff[i] += R[i][j] * Lambda[j] * dZij[j];
   }    
	
   flux.mass_flux       -= 0.5 * Diff[0];
   flux.momentum_flux.x -= 0.5 * Diff[1];
   flux.momentum_flux.y -= 0.5 * Diff[2];
   flux.momentum_flux.y -= 0.5 * Diff[3];
   flux.energy_flux     -= 0.5 * Diff[4]; 
   
   flux*=face_area;
   
   //cout<<flux.mass_flux<<" "<<flux.momentum_flux.x<<" "<<flux.momentum_flux.y<<" "<<flux.energy_flux<<endl;                                                    
}
