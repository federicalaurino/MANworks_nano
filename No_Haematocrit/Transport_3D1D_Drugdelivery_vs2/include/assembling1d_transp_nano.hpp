 /* -*- c++ -*- (enableMbars emacs c++ mode) */
/*==========================================================================
    "Computational models for nanoparticle transport in the vascular system"
    	       Master thesis in Mathematical Engineering
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
         Copyright (C) 2017 Annagiulia Tiozzo, Federica Laurino
============================================================================*/
/*! 
  @file   assembling1d_transp_nano.hpp
  @author Annagiulia Tiozzo <annagiulia92t@gmail.com>
	  Federica Laurino <federica.laurino@polimi.it>
  @date   April 2017.
  @brief  Miscelleanous assembly routines for the 1D nanoparticle transport problem.
 */
/** @defgroup asm Assembly routines */

#ifndef M3D1D_ASSEMBLING_1D_TRANSP_NANO_HPP_
#define M3D1D_ASSEMBLING_1D_TRANSP_NANO_HPP_
#include <defines.hpp>
#include <utilities.hpp>

namespace getfem {

/*! Build the adhesion term
    @f$ A=\int_{\Lambda}  M~c~b~d\sigma@f$
 */
/*!
	@param A         Computed mass matrix
	@param mim       The integration method to be used
	@param mf_c      The finite element method for the concentration @f$ c @f$
	@param mf_u  	 The finite element method for the adhesion coefficient
	@param adh	 The value for the adhesion coefficient
	@param rg        The region where to integrate

	@ingroup asm
 */

template<typename MAT, typename VEC>
void 
asm_network_nano_transp
	(MAT & A,  
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_u,
	 const VEC & adh,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_c.get_qdim() == 1 ,
		"invalid data mesh fem (Qdim=1 required)");
	
	generic_assembly 
	assem("m=data$1(#2);"
		  "t=comp(Base(#2).Base(#1).Base(#1));"
		  "M$1(#1,#1)+=t(k,:,:).m(k);"); 
	assem.push_mi(mim);
	assem.push_mf(mf_c);
	assem.push_mf(mf_u);
	assem.push_data(adh);
	assem.push_mat(A);
	assem.assembly(rg);

}

/*! Build the adhesion term
    @f$ A=\int_{\Lambda}  M~c~b~d\sigma@f$
 */
/*!
	@param A         Computed mass matrix
	@param mim       The integration method to be used
	@param mf_c      The finite element method for the concentration @f$ c @f$ and for the adhesion coefficient @f$ M @f$ 
	@param adh	 The value for the adhesion coefficient
	@param rg        The region where to integrate

	@ingroup asm
 */

template<typename MAT, typename VEC>
void 
asm_network_nano_transp
	(MAT & A,  
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const VEC & adh,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_c.get_qdim() == 1 ,
		"invalid data mesh fem (Qdim=1 required)");
	
	generic_assembly 
	assem("m=data$1(#1);"
		  "t=comp(Base(#1).Base(#1).Base(#1));"
		  "M$1(#1,#1)+=t(k,:,:).m(k);"); 
	assem.push_mi(mim);
	assem.push_mf(mf_c);
	assem.push_data(adh);
	assem.push_mat(A);
	assem.assembly(rg);

}


} /* end of namespace */

#endif
