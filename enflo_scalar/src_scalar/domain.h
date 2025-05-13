#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include <iostream>
#include <string>
#include "fparser.h"

//------------------------------------------------------------------------------
// Class to store limits of domain
//------------------------------------------------------------------------------
class Domain
{
   public:
      void add (std::string, std::string);
      void value (double &xmin, double &xmax, double &ymin, double &ymax, double &zmin, double &zmax);

   private:
      FParser Xmin;
      FParser Xmax;
      FParser Ymin;
      FParser Ymax;
      FParser Zmin;
      FParser Zmax;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void Domain::add (std::string variable, std::string fun)
{
   if(variable == "xmin")
      Xmin.FParse (fun);
   else if(variable == "xmax")
      Xmax.FParse (fun);
   else if(variable == "ymin")
      Ymin.FParse (fun);
   else if(variable == "ymax")
      Ymax.FParse (fun);      
   else if(variable == "zmin")
      Zmin.FParse (fun);
   else if(variable == "zmax")
      Zmax.FParse (fun); 
   else
   {
      std::cout << "Domain::add: Unknown domain limit " << variable << std::endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Evaluate domain limits
//------------------------------------------------------------------------------
inline
void Domain::value (double &xmin, double &xmax, double &ymin, double &ymax, double &zmin, double &zmax)
{
   double vals[] = {0, 0, 0, 0}; // DUMMY point. Not needed
   xmin = Xmin.Eval (vals);
   xmax = Xmax.Eval (vals);
   ymin = Ymin.Eval (vals);
   ymax = Ymax.Eval (vals);
   zmin = Zmin.Eval (vals);
   zmax = Zmax.Eval (vals);
}

#endif
