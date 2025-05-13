#ifndef __PRIMVAR_H__
#define __PRIMVAR_H__

#include "vec.h"

//------------------------------------------------------------------------------
//! Primary Variable class 
/*!
  This is describes primary variable as an object, along with various operators
  and functions
*/
//------------------------------------------------------------------------------
class PrimVar
{
   public:
      PrimVar(){density=0; pressure=0; velocity=0;};
      double density, pressure;
      Vector velocity;

      /*!
	  * Function saving the component-wise product of two primary variables into the 
	    current one
	  * @param[in] prim_var1
	  * @param[in] prim_var2
	  */
      PrimVar&  cdot(const PrimVar& prim_var1, const PrimVar& prim_var2); // componentwise multi
      
      /*!
	  * Function adding the component-wise product of two primary variables into the 
	    current one
	  * @param[in] prim_var1
	  * @param[in] prim_var2
	  */
      PrimVar&  addcdot(const PrimVar& prim_var1, const PrimVar& prim_var2); // componentwise multi
      
      /*!
	  * Function finding the component-wise product the absolute value of two primary variables 
	    and saving into the current one
	  * @param[in] prim_var1
	  * @param[in] prim_var2
	  */
      PrimVar&  adiff_eq(const PrimVar& prim_var1, const PrimVar& prim_var2); // componentwise multi
      
      /*!
	  * Operator multiplying the current primary variable by a scalar
	  * @param[in] scalar 
	  */
      PrimVar& operator*= (const double& scalar);
      
      /*!
	  * Operator assigning a scalar value to all primary variable components
	  * @param[in] scalar 
	  */
      PrimVar& operator=  (const double& scalar);
      
      /*!
	  * Operator adding a primary variable vector into the current one
	  * @param[in] prim_var
	  */
      PrimVar& operator+= (const PrimVar& prim_var);
      
      /*!
	  * Function setting the current primary variable as a linear sum of two 
	  * primitive variables
	  * @param[in] prim1
	  * @param[in] prim2
	  * @param[in] c1  
	  * @param[in] c2 
	  */
      PrimVar& equ(const PrimVar& prim1, const PrimVar& prim2, 
                   const double& c1, const double& c2);   
                   
      /*!
	  * Function setting the current primary variable as a linear sum of four 
	  * primitive variables
	  * @param[in] prim1
	  * @param[in] prim2
      * @param[in] prim3
	  * @param[in] prim4
	  * @param[in] c1  
	  * @param[in] c2 
      * @param[in] c3
	  * @param[in] c4
	  */
      PrimVar& equ(const PrimVar& prim1, const PrimVar& prim2, 
                   const PrimVar& prim3, const PrimVar& prim4,
                   const double& c1, const double& c2,
                   const double& c3, const double& c4);                       
      
      
      /*!
	  * Function adding to the current primary variable as a linear sum of four 
	  * primitive variables
	  * @param[in] prim1
	  * @param[in] prim2
      * @param[in] prim3
	  * @param[in] prim4
	  * @param[in] c1  
	  * @param[in] c2 
      * @param[in] c3
	  * @param[in] c4
	  */
      PrimVar& sadd(const PrimVar& prim1, const PrimVar& prim2, 
                   const PrimVar& prim3, const PrimVar& prim4,
                   const double& c1, const double& c2,
                   const double& c3, const double& c4);  
      
      /*!
	  * Function adding the squared (component-wise) primary variables into the given one 
	  * @param[in] prim
	  */ 
      PrimVar& sqadd(const PrimVar& prim);             
      
      /*!
	  * Function setting the current primary variable as the minimum (component-wise) 
	  * of between itself and another primary variable
	  * @param[in] p
	  */
      void min (const PrimVar& p);
      
      /*!
	  * Function setting the current primary variable as the maximum (component-wise) 
	  * of between itself and another primary variable
	  * @param[in] p
	  */
      void max (const PrimVar& p);
};


//------------------------------------------------------------------------------
// Multiply given primitive by scalar
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::operator*= (const double& scalar)
{
   density *= scalar;
   velocity    *= scalar; 
   pressure    *= scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Set a scalar value
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::operator= (const double& scalar)
{
   density = scalar;
   velocity    = scalar; 
   pressure    = scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Add primitive variable to given primitive variable
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::operator+= (const PrimVar& prim_var)
{
   density += prim_var.density;
   velocity    += prim_var.velocity;
   pressure    += prim_var.pressure;

   return *this;
}

//------------------------------------------------------------------------------
// Update PrimVar = min(PrimVar, p)
//------------------------------------------------------------------------------
inline
void PrimVar::min (const PrimVar& p)
{
   density = std::min(density, p.density);
   velocity.x  = std::min(velocity.x,  p.velocity.x);
   velocity.y  = std::min(velocity.y,  p.velocity.y);
   velocity.z  = std::min(velocity.z,  p.velocity.z);
   pressure    = std::min(pressure,    p.pressure);
}

//------------------------------------------------------------------------------
// Update PrimVar = max(PrimVar, p)
//------------------------------------------------------------------------------
inline
void PrimVar::max (const PrimVar& p)
{
   density = std::max(density, p.density);
   velocity.x  = std::max(velocity.x,  p.velocity.x);
   velocity.y  = std::max(velocity.y,  p.velocity.y);
   velocity.z  = std::max(velocity.z,  p.velocity.z);
   pressure    = std::max(pressure,    p.pressure);
}

//------------------------------------------------------------------------------
// Multiply two primitive variables componentwise
// Result is another primitive variable
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::cdot(const PrimVar& prim_var1, const PrimVar& prim_var2)
{

   density = prim_var1.density * prim_var2.density;
   velocity.x  = prim_var1.velocity.x  * prim_var2.velocity.x;
   velocity.y  = prim_var1.velocity.y  * prim_var2.velocity.y;
   velocity.z  = prim_var1.velocity.z  * prim_var2.velocity.z;
   pressure    = prim_var1.pressure    * prim_var2.pressure;

   return *this;
}


//------------------------------------------------------------------------------
// Multiply two primitive variables componentwise
// Result is added into another primitive variable
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::addcdot(const PrimVar& prim_var1, const PrimVar& prim_var2)
{

   density += prim_var1.density * prim_var2.density;
   velocity.x  += prim_var1.velocity.x  * prim_var2.velocity.x;
   velocity.y  += prim_var1.velocity.y  * prim_var2.velocity.y;
   velocity.z  += prim_var1.velocity.z  * prim_var2.velocity.z;
   pressure    += prim_var1.pressure    * prim_var2.pressure;

   return *this;
}

//------------------------------------------------------------------------------
// Find the absolute value of difference between two primitive variables componentwise
// and save into another primitive variable
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::adiff_eq(const PrimVar& prim_var1, const PrimVar& prim_var2)
{

   density     = std::fabs(prim_var1.density - prim_var2.density);
   velocity.x  = std::fabs(prim_var1.velocity.x  - prim_var2.velocity.x);
   velocity.y  = std::fabs(prim_var1.velocity.y  - prim_var2.velocity.y);
   velocity.z  = std::fabs(prim_var1.velocity.z  - prim_var2.velocity.z);
   pressure    = std::fabs(prim_var1.pressure    - prim_var2.pressure);

   return *this;
}

//------------------------------------------------------------------------------
// Add the square of difference between two primitive variables componentwise
// into another primitive variable
//------------------------------------------------------------------------------
// inline
// PrimVar& PrimVar::add_sqdiff(const PrimVar& prim_var1, const PrimVar& prim_var2)
// {
// 
//    density     += std::pow((prim_var1.density - prim_var2.density),2.0);
//    velocity.x  += std::pow((prim_var1.velocity.x  - prim_var2.velocity.x),2.0);
//    velocity.y  += std::pow((prim_var1.velocity.y  - prim_var2.velocity.y),2.0);
//    pressure    += std::pow((prim_var1.pressure    - prim_var2.pressure),2.0);
// 
//    return *this;
// }
// 
// //------------------------------------------------------------------------------
// // Find absolute value of difference between two primitive variables componentwise
// // and take the maximum value with current primitive variable
// //------------------------------------------------------------------------------
// inline
// PrimVar& PrimVar::max_adiff(const PrimVar& prim_var1, const PrimVar& prim_var2)
// {
// 
//    density     = std::max(std::fabs(prim_var1.density - prim_var2.density),density);
//    velocity.x  = std::max(std::fabs(prim_var1.velocity.x  - prim_var2.velocity.x),velocity.x);
//    velocity.y  = std::max(std::fabs(prim_var1.velocity.y  - prim_var2.velocity.y),velocity.y);
//    pressure    = std::max(std::fabs(prim_var1.pressure    - prim_var2.pressure),pressure);
// 
//    return *this;
// }


//------------------------------------------------------------------------------
// Adding two scaled primary variables and saving in given primary variable, 
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::equ(const PrimVar& prim1, const PrimVar& prim2,
                      const double& c1, const double& c2) 
{
   density  = prim1.density*c1  + prim2.density*c2;
   velocity.equ(prim1.velocity,prim2.velocity,c1,c2);
   pressure     = prim1.pressure*c1   + prim2.pressure*c2;
   
   return *this;
}

//------------------------------------------------------------------------------
// Adding four scaled primary variables and saving in given primary variable, 
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::equ(const PrimVar& prim1, const PrimVar& prim2, 
                      const PrimVar& prim3, const PrimVar& prim4,
                      const double& c1, const double& c2,
                      const double& c3, const double& c4) 
{
   density  = prim1.density*c1  + prim2.density*c2 + prim3.density*c3  + prim4.density*c4;
   velocity.equ(prim1.velocity,prim2.velocity,prim3.velocity,prim4.velocity,c1,c2,c3,c4);
   pressure  = prim1.pressure*c1   + prim2.pressure*c2 + prim3.pressure*c3   + prim4.pressure*c4;
   
   return *this;
}

//------------------------------------------------------------------------------
// Adding four scaled primary variables into given primary variable, 
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::sadd(const PrimVar& prim1, const PrimVar& prim2, 
                      const PrimVar& prim3, const PrimVar& prim4,
                      const double& c1, const double& c2,
                      const double& c3, const double& c4) 
{
   density  += prim1.density*c1  + prim2.density*c2 + prim3.density*c3  + prim4.density*c4;
   velocity.sadd(prim1.velocity,prim2.velocity,prim3.velocity,prim4.velocity,c1,c2,c3,c4);
   pressure  += prim1.pressure*c1   + prim2.pressure*c2 + prim3.pressure*c3   + prim4.pressure*c4;
   
   return *this;
}


//------------------------------------------------------------------------------
// Adding the squared (component-wise) primary variables into  given primary variable, 
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::sqadd(const PrimVar& prim)
{
   density  += prim.density*prim.density;
   velocity.sqadd(prim.velocity);
   pressure     += prim.pressure*prim.pressure;
   
   return *this;
}

#endif
