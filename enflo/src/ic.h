#ifndef __IC_H__
#define __IC_H__

#include <iostream>
#include <string>
#include <vector>
#include "fparser.h"
#include "ext_constants.h"
#include "vec.h"
#include "primvar.h"

//------------------------------------------------------------------------------
// Class to store initial condition functions
//------------------------------------------------------------------------------
class InitialCondition
{
   public:
      void add (std::string, std::string);
      void value (Vector& p, PrimVar & prim);

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
void InitialCondition::add (std::string variable, std::string fun)
{
   std::cout << fun <<std::endl;
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
      std::cout << "InitialCondition::add: Unknown variable " << variable << std::endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Evaluate primitive variables for given point
//------------------------------------------------------------------------------
inline
void InitialCondition::value (Vector& p, PrimVar &prim)
{

   double t = 0.0;
   double vals[] = {p.x, p.y, p.z, t};
   prim.density        = density.Eval (vals);
   prim.velocity.x     = xvelocity.Eval (vals);
   prim.velocity.y     = yvelocity.Eval (vals);
   prim.velocity.z     = zvelocity.Eval (vals);
   prim.pressure       = pressure.Eval (vals);
}

#endif
