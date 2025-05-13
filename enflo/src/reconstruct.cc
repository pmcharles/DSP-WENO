#include <cmath>
#include "fv.h"
#include <iomanip>
#include "limiter.h"

#define EPSILON  1.0e-14
#define KKK 1.0/3.0

using namespace std;

//------------------------------------------------------------------------------
// First order Reconstruct left and right states
//------------------------------------------------------------------------------
void FV::reconstruct_first (const std::vector< std::vector<double> > & state,
							double (& rec_statel)[NVAR],
							double (& rec_stater) [NVAR]) const
{
   int s = state.size();
   assert(s>=2);
   for(int i=0; i<NVAR; ++i)
   {
      rec_statel[i] = state[s/2-1][i];
      rec_stater[i] = state[s/2][i];
   }   
}



//------------------------------------------------------------------------------
// Second order Reconstruct left and right states (unlimited MUSCL)
//------------------------------------------------------------------------------
void FV::reconstruct_second_muscl (const std::vector< std::vector<double> > & state,
									 double (& rec_statel)[NVAR],
									 double (& rec_stater) [NVAR]) const
{
   double dstate, dstatel, dstater;
   int s = state.size();
   assert(s>=4);
   for(int i=0; i<NVAR; ++i)
   {
      dstate    = state[s/2][i] - state[s/2-1][i] ;
      dstatel   = state[s/2-1][i] - state[s/2-2][i] ;
      dstater   = state[s/2+1][i] - state[s/2][i] ;
      
      // left state
      rec_statel[i] = state[s/2-1][i] + (dstatel*(1.0-KKK) + dstate*(1.0+KKK))*0.25;
      
      // right state
      rec_stater[i] = state[s/2][i] - (dstater*(1.0-KKK) + dstate*(1.0+KKK))*0.25;
   }
}

//------------------------------------------------------------------------------
// Reconstruct left and right states (limited MUSCL)
//------------------------------------------------------------------------------
void FV::reconstruct_muscl_limited(const ReconstructionScheme reconstruct_scheme,
       							   const double& h,
								   const std::vector< std::vector<double> > & state,
								   double (& rec_statel)[NVAR],
								   double (& rec_stater) [NVAR]) const
{
   double dstate, dstatel, dstater;
   double sl,sr;
   int s = state.size();
   assert(s>=4);
   for(int i=0; i<NVAR; ++i)
   {
      dstate    = state[s/2][i] - state[s/2-1][i] ;
      dstatel   = state[s/2-1][i] - state[s/2-2][i] ;
      dstater   = state[s/2+1][i] - state[s/2][i] ;
      
      // left state
      sl = limited_slope(dstatel,dstate,reconstruct_scheme,h);
      rec_statel[i] = state[s/2-1][i] + sl;
      
      // right state
      sr = limited_slope(dstater,dstate,reconstruct_scheme,h);
      rec_stater[i] = state[s/2][i] - sr;
   }
}

//------------------------------------------------------------------------------
// Reconstruct left and right states using point wise ENO
//------------------------------------------------------------------------------
void FV::reconstruct_eno_pw(const int eno_level,
							const double& h,
							const std::vector< std::vector<double> > & state,
							double (& rec_statel)[NVAR],
							double (& rec_stater) [NVAR]) const
{
   int s = state.size();
   assert(s>=2*eno_level);
   double state_comp[2*eno_level-1];
   for(int i=0; i<NVAR; ++i)
   {
      //left state
      for(int j=0; j<2*eno_level-1; ++j)
      {
         state_comp[j] = state[s/2 - eno_level + j][i];
      }
      rec_statel[i]=eno_pw(eno_level,state_comp);
      
      // right state
      for(int j=0; j<2*eno_level-1; ++j)
      {
         state_comp[2*eno_level-2 -j] = state[s/2 - eno_level + 1 + j][i];
      }
      rec_stater[i]=eno_pw(eno_level,state_comp);
   }
}

//------------------------------------------------------------------------------
// Reconstruct left and right states using point wise SP-WENO
//------------------------------------------------------------------------------
void FV::reconstruct_sp_weno(const double& h,
							const std::vector< std::vector<double> > & state,
							double (& rec_statel)[NVAR],
							double (& rec_stater) [NVAR]) const
{
   int s = state.size();
   assert(s>=4);
   double state_comp[4];
   for(int i=0; i<NVAR; ++i)
   {
      // All four cell values for component i
      for(int j=0; j<4; ++j)
      {
         state_comp[j] = state[s/2 - 2 + j][i];
      }

      sp_weno_rec(state_comp,rec_statel[i],rec_stater[i],param.spweno_corr_type);
   }
}

//------------------------------------------------------------------------------
// Reconstruct left and right states using point wise SP-WENO
//------------------------------------------------------------------------------
void FV::reconstruct_sp_weno_dl(const std::vector< std::vector<double> > & state,
							double (& rec_statel)[NVAR],
							double (& rec_stater) [NVAR]) const
{
   int s = state.size();
   assert(s>=4);
   double state_comp[4];
   for(int i=0; i<NVAR; ++i)
   {
      // All four cell values for component i
      for(int j=0; j<4; ++j)
      {
         state_comp[j] = state[s/2 - 2 + j][i];
      }
      if (abs(state_comp[2] - state_comp[1]) < 1e-15)
         rec_statel[i] = 0.0, rec_stater[i] = 0.0;
      else
         net.forward(state_comp,rec_statel[i],rec_stater[i]);
   }
}

//------------------------------------------------------------------------------
// Computed limited slope 
//------------------------------------------------------------------------------
double FV::limited_slope (const double& ul, 
						  const double& ur,
						  const ReconstructionScheme reconstruct_scheme,
						  const double& h) const
{
   double result;
   switch(reconstruct_scheme)
   {
      case tvd_minmod:
         result = 0.5*minmod(ul,ur);
         break;
      
      case minabs_lim:
         result = 0.5*minabs(ul,ur);
         break;   
         
      default:
         cout << "reconstruct: unknown reconstruction scheme = " 
              << reconstruct_scheme << endl;
         exit (0);
   }
   

   return result;
}

//------------------------------------------------------------------------------
// Reconstruct left and right states 
// V.size() == 2*ng_cell
//------------------------------------------------------------------------------
void FV::reconstruct (const ReconstructionScheme reconstruct_scheme,
                      const double h,
					  const std::vector< std::vector<double> > & state,
					  double (& rec_statel)[NVAR],
					  double (& rec_stater) [NVAR]) const
{
   switch(reconstruct_scheme)
   {
      // First order
      case first:
         reconstruct_first (state,rec_statel,rec_stater);
         break;

      // Second order
      case second:
         reconstruct_second_muscl (state,rec_statel,rec_stater);
         break;
         
      // MUSCL with minmod limiter
      case tvd_minmod:
         reconstruct_muscl_limited (reconstruct_scheme,h,state,rec_statel,rec_stater);
         break;
         
      // MUSCL with minabs limiter
      case minabs_lim:
         reconstruct_muscl_limited (reconstruct_scheme,h,state,rec_statel,rec_stater);
         break;   
         
      // ENO
      case eno:
         reconstruct_eno_pw(param.eno_level,h,state,rec_statel,rec_stater);
         break;
         
      // SP-WENO
      case sp_weno:
         reconstruct_sp_weno(h,state,rec_statel,rec_stater);
         break;      
      
      // SP-WENO-DL
      case sp_weno_dl:
		   reconstruct_sp_weno_dl(state,rec_statel,rec_stater);
		   break;
     
      default:
         cout<<"reconstruct: unknown reconstruction scheme = "<<reconstruct_scheme<<endl;
         exit(0);
   }
}
