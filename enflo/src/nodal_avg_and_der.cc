#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Finds primitive average at node and derivatives of velocity
// Using multi-dimensional operators
// (i,j,k) --> corresponds to node (i+h/2,j+h/2,k+h/2) and will 
// represent the index for the primitive variable at cell-center 
// [ng_cell + i][ng_cell + j][ng_cell +k]
// Algorithm needs one ghost layer around mesh
//------------------------------------------------------------------------------

void FV::nodal_avg_and_vel_derivatives (int i, int j, int k,
                                        PrimVar &prim_avg, Vector &Txyz_vel,
                                        Vector &Tyxz_vel, Vector &Tzxy_vel) const
{  
   assert(i>=-1); assert(i<Nx_cell);
   assert(j>=-1); assert(j<Ny_cell);
   assert(k>=-1); assert(k<Nz_cell);
           
   double fact1 = 1.0/8.0;        
   
   prim_avg.equ(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1],
                primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1], 
                primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1],
                primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k   ],
                fact1, fact1, fact1, fact1);
   prim_avg.sadd(primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k +1],
                 primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k   ],
                 primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k   ],
                 primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k   ],
                 fact1, fact1, fact1, fact1);    
   
   double factx = 1.0/(4.0*dx), facty = 1.0/(4.0*dy), factz = 1.0/(4.0*dz);
   

   Txyz_vel.equ(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1].velocity,
                primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k   ].velocity, 
                primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1].velocity,
                primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k   ].velocity,
                factx, factx, factx, factx);
   Txyz_vel.sadd(primitive[ng_cell + i][ng_cell + j +1][ng_cell +k +1].velocity,
                 primitive[ng_cell + i][ng_cell + j +1][ng_cell +k   ].velocity, 
                 primitive[ng_cell + i][ng_cell + j   ][ng_cell +k +1].velocity,
                 primitive[ng_cell + i][ng_cell + j   ][ng_cell +k   ].velocity,
                 -factx, -factx, -factx, -factx);
        
              
   Tyxz_vel.equ(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1].velocity,
                primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k   ].velocity, 
                primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1].velocity,
                primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k   ].velocity,
                facty, facty, facty, facty);
   Tyxz_vel.sadd(primitive[ng_cell + i +1][ng_cell + j][ng_cell +k +1].velocity,
                 primitive[ng_cell + i +1][ng_cell + j][ng_cell +k   ].velocity, 
                 primitive[ng_cell + i   ][ng_cell + j][ng_cell +k +1].velocity,
                 primitive[ng_cell + i   ][ng_cell + j][ng_cell +k   ].velocity,
                 -facty, -facty, -facty, -facty);
                      
   Tzxy_vel.equ(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1].velocity,
                primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1].velocity, 
                primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1].velocity,
                primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k +1].velocity,
                factz, factz, factz, factz);
   //cout<<primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1].velocity.x<<" "
       //<<primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k +1].velocity.x<<" "
       //<<primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k +1].velocity.x<<" "
       //<<primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k +1].velocity.x<<" "
       //<<primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k +1].velocity.x<<" "
       //<<primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k].velocity.x<<" "
       //<<primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k].velocity.x<<" "
       //<<primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k].velocity.x<<" ";
        
   Tzxy_vel.sadd(primitive[ng_cell + i +1][ng_cell + j +1][ng_cell +k].velocity,
                 primitive[ng_cell + i +1][ng_cell + j   ][ng_cell +k].velocity, 
                 primitive[ng_cell + i   ][ng_cell + j +1][ng_cell +k].velocity,
                 primitive[ng_cell + i   ][ng_cell + j   ][ng_cell +k].velocity,
                 -factz, -factz, -factz, -factz);
}
