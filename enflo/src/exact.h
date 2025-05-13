#ifndef __EXACT_H__
#define __EXACT_H__

#include <iostream>
#include <string>
#include <vector>
#include "fparser.h"
#include "ext_constants.h"
#include "vec.h"
#include "primvar.h"

//------------------------------------------------------------------------------
// Class to store Exact solution functions
//------------------------------------------------------------------------------
class Exact
{
   public:
      void    add (std::string, std::string);
      void value (Vector& p,double &t, PrimVar & prim);

   private:
      FParser density;
      FParser xvelocity;
      FParser yvelocity;
      FParser zvelocity;
      FParser pressure;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void Exact::add (std::string variable, std::string fun)
{
   if(variable == "density")
      density.FParse (fun);
   else if(variable == "xvelocity")
      xvelocity.FParse (fun);
   else if(variable == "yvelocity")
      yvelocity.FParse (fun);  
   else if(variable == "zvelocity")
      zvelocity.FParse (fun);      
   else if(variable == "pressure")
      pressure.FParse (fun);
   else
   {
      std::cout << "Exact::add: Unknown variable " << variable << std::endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Evaluate primitive variables for given point
//------------------------------------------------------------------------------
inline
void Exact::value (Vector& p,double &t, PrimVar &prim)
{

   double vals[] = {p.x, p.y, p.z, t};
   prim.density        = density.Eval (vals);
   prim.velocity.x     = xvelocity.Eval (vals);
   prim.velocity.y     = yvelocity.Eval (vals);
   prim.velocity.z     = zvelocity.Eval (vals);
   prim.pressure       = pressure.Eval (vals);
}

#endif
