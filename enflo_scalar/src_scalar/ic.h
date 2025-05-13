#ifndef __IC_H__
#define __IC_H__

#include <iostream>
#include <string>
#include "fparser.h"
#include "ext_constants.h"
#include "vec.h"

//------------------------------------------------------------------------------
// Class to store initial condition functions
//------------------------------------------------------------------------------
class InitialCondition
{
   public:
      void add (std::string);
      void value (Vector& p,double & soln);

   private:
      FParser initial;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void InitialCondition::add (std::string fun)
{
   initial.FParse (fun);
}

//------------------------------------------------------------------------------
// Evaluate solution at given point
//------------------------------------------------------------------------------
inline
void InitialCondition::value (Vector& p, double &soln)
{
   double t = 0.0;
   double vals[] = {p.x, p.y, p.z, t};
   soln        = initial.Eval (vals);
}

#endif
