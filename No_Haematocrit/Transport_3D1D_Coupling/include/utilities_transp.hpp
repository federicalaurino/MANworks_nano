/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   utilities_transp.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   November 2016.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */
#ifndef M3D1D_UTILITIES_TRANSP_HPP_
#define M3D1D_UTILITIES_TRANSP_HPP_

#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_mesh.h>
#include <gmm/gmm.h>

namespace gmm {

//! Componentwise product of two vectors:
//!  B[i]= A[i]* B[i]

template
<typename VEC>
void scale (const VEC & A, VEC & B){

//calculate the length of the two vectors

int lengthA = gmm::mat_nrows(gmm::col_vector(A));
int lengthB = gmm::mat_nrows(gmm::col_vector(B));
GMM_ASSERT1(lengthA==lengthB, "impossible to scale the vectors: different lengths");

for(int i=0; i<lengthA; i++)
  B[i]= A[i]*B[i];
  
} 


} /* end of namespace */

#endif
