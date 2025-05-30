#ifndef __ENTVAR_H__
#define __ENTVAR_H__

#include "vec.h"

//------------------------------------------------------------------------------
//! Entropy Variable class 
/*!
  This is describes entropy variable as an object, along with various operators
  and functions
*/
//------------------------------------------------------------------------------
class EntVar
{
   public:
      EntVar(){entf = 0; entl=0; entv=0;};
      double entf, /*!< (f)irst component of the entropy variable corresponding to
                        density */ 
             entl; /*!< (l)ast component of the entropy variable corresponding to energy */
      Vector entv; /*!< components of the entropy variable corresponding to the momentum 
                        vector  */
      
      /*!
	  * Operator multiplying the current entropy variable by a scalar
	  * @param[in] scalar 
	  */
      EntVar& operator*= (const double& scalar);
      
      /*!
	  * Operator assigning a scalar value to all entropy variable components
	  * @param[in] scalar 
	  */
      EntVar& operator=  (const double& scalar);
      
      /*!
	  * Operator adding an entropy variable vector into the current one
	  * @param[in] ent_var
	  */
      EntVar& operator+= (const EntVar& ent_var);
      
      /*!
	  * Function setting the current entropy variable as a linear sum of a two 
	  * entropy variables
	  * @param[in] ent1
	  * @param[in] ent2
	  * @param[in] c1  
	  * @param[in] c2 
	  */
      EntVar& equ(const EntVar& ent1, const EntVar& ent2, 
                   const double& c1, const double& c2);
                   
      /*!
	  * Function setting the current entropy variable as a linear sum of a three 
	  * entropy variables
	  * @param[in] ent1
	  * @param[in] ent2
	  * @param[in] ent3
	  * @param[in] c1  
	  * @param[in] c2 
	  * @param[in] c3
	  */             
      EntVar& equ(const EntVar& ent1, const EntVar& ent2, const EntVar& ent3, 
                   const double& c1, const double& c2, const double& c3);     
                   
      /*!
	  * Function setting the current entropy variable as a linear sum of four 
	  * entropy variables
	  * @param[in] ent1
	  * @param[in] ent2
      * @param[in] ent3
	  * @param[in] ent4
	  * @param[in] c1  
	  * @param[in] c2 
      * @param[in] c3
	  * @param[in] c4
	  */
      EntVar& equ(const EntVar& ent1, const EntVar& ent2, 
                   const EntVar& ent3, const EntVar& ent4,
                   const double& c1, const double& c2,
                   const double& c3, const double& c4);                       
      
      
      /*!
	  * Function adding to the current entropy variable as a linear sum of four 
	  * entropy variables
	  * @param[in] ent1
	  * @param[in] ent2
      * @param[in] ent3
	  * @param[in] ent4
	  * @param[in] c1  
	  * @param[in] c2 
      * @param[in] c3
	  * @param[in] c4
	  */
      EntVar& sadd(const EntVar& ent1, const EntVar& ent2, 
                   const EntVar& ent3, const EntVar& ent4,
                   const double& c1, const double& c2,
                   const double& c3, const double& c4);                      
      
      /*!
	  * Function setting the current entropy variable as the minimum (component-wise) 
	  * of between itself and another entropy variable
	  * @param[in] p
	  */
      void min (const EntVar& p);
      
      /*!
	  * Function setting the current entropy variable as the maximum (component-wise) 
	  * of between itself and another entropy variable
	  * @param[in] p
	  */
      void max (const EntVar& p);
      
      /*!
	  * Function setting the current entropy variable as the maximum (component-wise) 
	  * of between itself and another entropy variable in absolute values
	  * @param[in] p
	  */
      void maxabs (const EntVar& p);
};

//------------------------------------------------------------------------------
// Multiply given Entropy by scalar
//------------------------------------------------------------------------------
inline
EntVar& EntVar::operator*= (const double& scalar)
{
   entf *= scalar;
   entv    *= scalar; 
   entl    *= scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Set a scalar value
//------------------------------------------------------------------------------
inline
EntVar& EntVar::operator= (const double& scalar)
{
   entf = scalar;
   entv    = scalar; 
   entl    = scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Add Entropy variable to given Entropy variable
//------------------------------------------------------------------------------
inline
EntVar& EntVar::operator+= (const EntVar& ent_var)
{
   entf += ent_var.entf;
   entv    += ent_var.entv;
   entl    += ent_var.entl;

   return *this;
}

//------------------------------------------------------------------------------
// Update EntVar = min(EntVar, p)
//------------------------------------------------------------------------------
inline
void EntVar::min (const EntVar& p)
{
   entf = std::min(entf, p.entf);
   entv.x  = std::min(entv.x,  p.entv.x);
   entv.y  = std::min(entv.y,  p.entv.y);
   entv.z  = std::min(entv.z,  p.entv.z);
   entl    = std::min(entl,    p.entl);
}

//------------------------------------------------------------------------------
// Update EntVar = max(EntVar, p)
//------------------------------------------------------------------------------
inline
void EntVar::max (const EntVar& p)
{
   entf = std::max(entf, p.entf);
   entv.x  = std::max(entv.x,  p.entv.x);
   entv.y  = std::max(entv.y,  p.entv.y);
   entv.z  = std::max(entv.z,  p.entv.z);
   entl    = std::max(entl,    p.entl);
}

//------------------------------------------------------------------------------
// Update EntVar = max(|EntVar|, |p|)
//------------------------------------------------------------------------------
inline
void EntVar::maxabs (const EntVar& p)
{
   entf    = std::max(std::fabs(entf), std::fabs(p.entf));
   entv.x  = std::max(std::fabs(entv.x),  std::fabs(p.entv.x));
   entv.y  = std::max(std::fabs(entv.y),  std::fabs(p.entv.y));
   entv.z  = std::max(std::fabs(entv.z),  std::fabs(p.entv.z));
   entl    = std::max(std::fabs(entl),    std::fabs(p.entl));
}


//------------------------------------------------------------------------------
// Adding 2 scaled entropy variables into  given entropy variable, 
//------------------------------------------------------------------------------
inline
EntVar& EntVar::equ(const EntVar& ent1, const EntVar&  ent2, 
                     const double& c1, const double& c2)
{
   entf  = ent1.entf*c1  + ent2.entf*c2;
   entv.equ(ent1.entv,ent2.entv,c1,c2);
   entl  = ent1.entl*c1   + ent2.entl*c2;
   
   return *this;
}

//------------------------------------------------------------------------------
// Adding 3 scaled entropy variables into  given entropy variable, 
//------------------------------------------------------------------------------
inline
EntVar& EntVar::equ(const EntVar& ent1, const EntVar&  ent2, const EntVar&  ent3,
                     const double& c1, const double& c2, const double& c3)
{
   entf  = ent1.entf*c1  + ent2.entf*c2 + ent3.entf*c3;
   entv.equ(ent1.entv,ent2.entv,ent3.entv,c1,c2,c3);
   entl  = ent1.entl*c1 + ent2.entl*c2 + ent3.entl*c3;
   
   return *this;
}

//------------------------------------------------------------------------------
// Adding four scaled entropy variables and saving in given entropy variable, 
//------------------------------------------------------------------------------
inline
EntVar& EntVar::equ(const EntVar& ent1, const EntVar& ent2, 
                      const EntVar& ent3, const EntVar& ent4,
                      const double& c1, const double& c2,
                      const double& c3, const double& c4) 
{
   entf  = ent1.entf*c1  + ent2.entf*c2 + ent3.entf*c3  + ent4.entf*c4;
   entv.equ(ent1.entv,ent2.entv,ent3.entv,ent4.entv,c1,c2,c3,c4);
   entl  = ent1.entl*c1   + ent2.entl*c2 + ent3.entl*c3   + ent4.entl*c4;
   
   return *this;
}

//------------------------------------------------------------------------------
// Adding four scaled entropy variables into given entropy variable, 
//------------------------------------------------------------------------------
inline
EntVar& EntVar::sadd(const EntVar& ent1, const EntVar& ent2, 
                      const EntVar& ent3, const EntVar& ent4,
                      const double& c1, const double& c2,
                      const double& c3, const double& c4)
{
   entf  += ent1.entf*c1  + ent2.entf*c2 + ent3.entf*c3  + ent4.entf*c4;
   entv.sadd(ent1.entv,ent2.entv,ent3.entv,ent4.entv,c1,c2,c3,c4);
   entl  += ent1.entl*c1   + ent2.entl*c2 + ent3.entl*c3   + ent4.entl*c4;
   
   return *this;
}



#endif
