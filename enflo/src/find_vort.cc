#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Computing vortcity at nodes
//------------------------------------------------------------------------------
void FV::find_vorticity()
{
   PrimVar prim_avg;
   Vector Txyz_vel, Tyxz_vel, Tzxy_vel;
   
   // Simulataneous loop for cell-centers and nodes
   // (i,j,k) --> (i,j,k) cell-center 
   // (i,j,k) --> (i+h/1,j+h/2k+h/2)
   
   for(int i=-1; i<Nx_cell; ++i)
      for(int j=-1; j<Ny_cell; ++j)
         for(int k=-1; k<Nz_cell; ++k)
		 {
            nodal_avg_and_vel_derivatives (i, j, k, prim_avg, Txyz_vel, Tyxz_vel,Tzxy_vel);
            //cout<<Tzxy_vel.y<<" "<<Tzxy_vel.x;
            vorticity[i+1][j+1][k+1].x = Tyxz_vel.z - Tzxy_vel.y;
            vorticity[i+1][j+1][k+1].y = Tzxy_vel.x - Txyz_vel.z;
            vorticity[i+1][j+1][k+1].z = Txyz_vel.y - Tyxz_vel.x;
            
		 }
}
