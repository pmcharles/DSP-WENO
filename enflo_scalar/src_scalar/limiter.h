#ifndef __LIMITER_H__
#define __LIMITER_H__

#include <cmath>

#define SIGN(a) (((a)<0) ? -1:1)

const double ENO2_COEF[][2] = {{-0.5,1.5},
							   {0.5,0.5}};
const double ENO3_COEF[][3] = {{0.375,-1.25,1.875},
							   {-0.125,0.75,0.375},
							   {0.375,0.75,-0.125}}; 

const double eps = 1.0e-13;							                    

using namespace std;

//------------------------------------------------------------------------------
// Minmod of two number
//------------------------------------------------------------------------------
inline
double minmod (const double& a,
               const double& b)
{
   double result;
   
   if (a*b > 0.0)
   {
      result  = min(fabs(a), fabs(b));
      result *= SIGN(a);
   }
   else
      result = 0.0;
   
   return result;
   
}

//------------------------------------------------------------------------------
// MINABS
//------------------------------------------------------------------------------
inline
double minabs (const double& a,
               const double& b)
{
   double result;
   
   if (fabs(a) > fabs(b))
   {
      result  = b;
   }
   else
      result = a;
   
   return result;
   
}

//------------------------------------------------------------------------------
// Pointwise ENO reconstruction
//------------------------------------------------------------------------------
inline
double eno_pw (const int &eno_level, const double state[])
{
   double result = 0;
   
   int r = 0; // shift in stencil position
   double dd[eno_level][2*eno_level-1];
   
   
   //Calculating divided differences
   for(int j=0; j<2*eno_level-1; ++j)
      dd[0][j] = state[j];
   
   for(int i=1; i<eno_level; ++i)
      for(int j=0; j<2*eno_level-1- i; ++j)
         dd[i][j] = dd[i-1][j+1] - dd[i-1][j];
   
   // Finding r shift
   int s = eno_level-1;
   for(int i=1; i<eno_level; ++i)
   {
      if(abs(dd[i][s-r]) > abs(dd[i][s-r-1]))
        r++;
   }
   
   // Interpolating
   for(int i=0; i<eno_level; ++i)
   {
      if(eno_level == 2)
         result += ENO2_COEF[s-r][i]*state[s-r+i];
      else if(eno_level == 3)   
         result += ENO3_COEF[s-r][i]*state[s-r+i];
   }
   return result;
   
}

//------------------------------------------------------------------------------
// Needed for sp-weno
//------------------------------------------------------------------------------
inline
double fp (const double &theta0, const double &theta1)
{
   double psi, val;
   
   if(fabs(theta0 - 1.0)<eps)
     psi = 0.0;
   else
     psi = (1-theta1)/(1-theta0);
     
   if(fabs(theta0 - 1.0)>eps && fabs(psi + 1.0)>eps)
      val = 1.0/(1.0 + psi);
   else
      val = 1.0;    
   
   return val;       
}

// inline
// double Cp (const double &theta0, const double &theta1)
// {
//    double psi, val = 1.0/8.0; 
//    
//    if(theta0 == 1.0)//(fabs(theta0 - 1.0)<eps)
//      psi = 0.0;
//    else
//      psi = (1-theta1)/(1-theta0);
//      
//    if(theta0!=1.0 && psi < 0 && psi !=-1)//(fabs(psi + 1.0)>eps && (theta0-1.0)>eps && psi < 0)
//       if((theta0<-1 && theta1 > 1) || (theta1<-1 && theta0 > 1))
//          val = 1.0/8.0;
//       else   
//          val = (fp(theta0,theta1)/(pow(fp(theta0,theta1),2.0) + pow(fp(theta1,theta0),2.0)))/8.0;
//    else if(theta0!=1.0 && psi == -1.0)//((theta0-1.0)>eps && fabs(psi + 1.0)<eps)
//       val = 0.0;  
//    else if(theta0 == 1.0 ||(psi >=0 && fabs(theta0)<=1.0))//((theta0-1.0)<eps || (psi > -eps && fabs(theta0)<1+eps))
//       val = -3.0/8.0;
//    else if(psi >=0 && fabs(theta0)>1.0) //(psi > -eps && fabs(theta0)>1)
//       val = 1.0/8.0; 
//    else
//    {
//       cout<<"WARNING: UNKNOWN SP-WENO REGION ENCOUNTERED!!"<<endl;  
// //       cout<<theta0<<" "<<theta1<<endl;
// //       cout<<del0<<" "<<del1<<" "<<del2<<endl;
// //       cout<<state[0]<<" "<<state[1]<<" "<<state[2]<<" "<<state[3]<<endl;
//    }           
//    
//    return val;       
// }

inline
void Cp (const double &theta0, const double &theta1,const double &del1,const double &scale, double &originalC, double &correction)
{
   double psi; 
   
   if(fabs(theta0 - 1.0)<eps || fabs(theta1 - 1.0)<eps)
   psi = 0.0;
   else
     psi = (1-theta1)/(1-theta0);
     
   correction = 0.0;  
   if(fabs(psi + 1.0)>eps && fabs(theta0-1.0)>eps && psi < 0)
   {
      originalC = (fp(theta0,theta1)/(pow(fp(theta0,theta1),2.0) + pow(fp(theta1,theta0),2.0)))/8.0;
      correction += - 0.25*pow(min(fabs(del1/scale),fabs(del1)),3.0)/(1.0-theta0); 
   }   
   else if(fabs(theta0-1.0)>eps && fabs(psi + 1.0)<eps)
   {
      originalC = 0.0;
   }     
   else if(fabs(theta0-1.0)<eps || (psi > -eps && fabs(theta0)<1+eps))
   {
      originalC = -3.0/8.0;
   }   
   else if(psi > -eps && fabs(theta0)>1)
   {
      originalC = 1.0/8.0; 
   }
   else
   {
      cout<<"WARNING: UNKNOWN SP-WENO REGION ENCOUNTERED!!"<<endl;  
      cout<<theta0<<" "<<theta1<<endl;
      // cout<<del0<<" "<<del1<<" "<<del2<<endl;
//       cout<<state[0]<<" "<<state[1]<<" "<<state[2]<<" "<<state[3]<<endl;
   }                  
}

//------------------------------------------------------------------------------
// Pointwise WENO reconstruction
//------------------------------------------------------------------------------
inline
void sp_weno_rec (const double state[], double &recl, double &recr, const int &corr_type)
{
   double theta_i_p, theta_ip1_m, scale;
   double del0 = state[1] - state[0];
   double del1 = state[2] - state[1];
   double del2 = state[3] - state[2];
   
   double w0, w0t;
   
   if(fabs(del1)<eps)//(fabs(del1) < eps)
   {
      recl = state[1];
      recr = state[2];
   }
   else
   {
      theta_i_p = del0/del1;
      theta_ip1_m = del2/del1;
      scale = 0.5*(fabs(state[1]) + fabs(state[2]));
      
      double C1=0,C2=0,C1_corr=0,C2_corr=0;
      
      //cout<<theta_i_p<<" "<<theta_ip1_m<<endl;
      Cp(theta_i_p, theta_ip1_m,del1,scale,C1,C1_corr);
      Cp(theta_ip1_m,theta_i_p,del1,scale,C2,C2_corr);
      
      if(corr_type > 0)
      {
          double C1n = C1 + C1_corr;
          double C2n = C2 + C2_corr;
          // THIS IS SP-WENOc
          if(corr_type == 1)
          {
              C1 = C1n;
              C2 = C2n;
          }
          // THIS IS EXPERIMENTAL MODIFICATION OF SP-WENOc
          // AS OF NOW IT IS THE SAME AS ORIGINAL SP-WENO
          else if(corr_type == 2)
          {
             //if((C1n < -3.0/8.0 && C2n > 1.0/8.0) ||(C2n < -3.0/8.0 && C1n > 1.0/8.0) &&
             // !((theta_i_p < -1.0 && theta_ip1_m >1.0) || (theta_ip1_m < -1.0 && theta_i_p >1.0)))
              //{
                  //C1 = C1n;
                  //C2 = C2n;
              //}
          }
          C1 = min(max(C1,-3.0/8.0),1.0/8.0);
          C2 = min(max(C2,-3.0/8.0),1.0/8.0);
      }
      
      w0 = 3.0/4.0 + 2*C1;
      w0t = 1.0/4.0 - 2*C2;
        
      recl = w0*0.5*(state[1] + state[2]) + (1-w0)*(-0.5*state[0] + 1.5*state[1]);
      recr = (1-w0t)*0.5*(state[1] + state[2]) + w0t*(-0.5*state[3] + 1.5*state[2]);
   }
   
   
}


#endif
