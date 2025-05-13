#ifndef __MISC_FUNC_H__
#define __MISC_FUNC_H__

#include <cmath>

#define SIGN(a) (((a)<0) ? -1:1)

using namespace std;


//------------------------------------------------------------------------------
// log average
//------------------------------------------------------------------------------
inline
double logavg(double a, double b)
{
   double xi = b/a;
   double f = (xi - 1.0) / (xi + 1.0);
   double u = f * f;

   double F;
   if (u < 1.0e-2)
   {
      double u2 = u * u;
      double u3 = u2 * u;
      F = 1.0 + u/3.0 + u2/5.0 + u3/7.0;
   }
   else
      F = log(xi)/2.0/f;

   return 0.5*(a+b)/F;
}

#endif
