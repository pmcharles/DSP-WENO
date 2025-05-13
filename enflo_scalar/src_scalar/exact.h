#ifndef __EXACT_H__
#define __EXACT_H__

#include <iostream>
#include <string>
#include "fparser.h"
#include "ext_constants.h"
#include "vec.h"

//------------------------------------------------------------------------------
// Class to store Exact solution functions
//------------------------------------------------------------------------------
class Exact
{
   public:
      void    add (std::string);
      void value (Vector& p,double &t, double & soln);

   private:
      FParser exact;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void Exact::add (std::string fun)
{
   exact.FParse (fun);
}

//------------------------------------------------------------------------------
// Evaluate exact soln at a  given point
//------------------------------------------------------------------------------
inline
void Exact::value (Vector& p,double &t, double &soln)
{

   double vals[] = {p.x, p.y, p.z, t};
   soln        = exact.Eval (vals);
}

#endif
