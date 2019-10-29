 /* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   assembling1d.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Miscelleanous assembly routines for the 1D transport problem.
 */
/** @defgroup asm Assembly routines */

#ifndef M3D1D_ASSEMBLING_1D_TRANSP_HPP_
#define M3D1D_ASSEMBLING_1D_TRANSP_HPP_
#include <defines.hpp>
#include <node.hpp>
#include <utilities.hpp>
#include <utilities_transp_nano.hpp>

namespace getfem {

//! Build the mass and diffusion matrices for the 1D Poiseuille's problem,
//! @f$ T = \int_{\Lambda} ~c~b~ds @f$ and
//! @f$ D = \int_{\Lambda} \nabla c \cdot \nabla b~ds @f$
/*!
	@param T         Computed mass matrix
	@param D         Computed diffusion matrix
	@param mim       The integration method to be used
	@param mf_c      The finite element method for the concentration @f$ c @f$
	@param mf_data   The finite element method for the diffusion coefficient
	@param diff  	 The diffusion coefficient for D
	@param rg        The region where to integrate

	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_network_poiseuille_transp
	(MAT & D, MAT & T, 
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data,
	 VEC & diff,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_c.get_qdim() == 1 ,
		"invalid data mesh fem (Qdim=1 required)");
	//build mass matrix Tv for time derivative
	getfem::asm_mass_matrix(T, mim, mf_c, rg);
	masslumping(T);
	// Build the diffusion matrix Dv
	getfem::asm_stiffness_matrix_for_laplacian(D,mim,mf_c,mf_data, diff, rg);

}

//! Build the advection matrix for the 1D Poiseuille's problem,
//! @f$ B = \int_{\Lambda} ~U \cdot \nabla c~b~ds @f$ and
/*!
	@param D         Computed advection matrix
	@param mim       The integration method to be used
	@param mf_c      The finite element method for the concentration @f$ c @f$
	@param mf_data   The finite element method for the tangent vector @f$ \Lambda @f$
	@param mf_u      The finite element method for the velocity vector @f$ U @f$
	@param U      	 The velocity vector field
	@param lambdax   First cartesian component of the tangent versor  @f$ \mathbf{\lambda} @f$
	@param lambday   Second cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param lambdaz   Third cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param rg        The region where to integrate

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void 
asm_advection_network
	(MAT & D,
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data,
	 const mesh_fem & mf_u,
	 const VEC & U,
	 const VEC & lambdax, const VEC & lambday, const VEC & lambdaz,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 	
	 
	 {

generic_assembly 
	assem("l1=data$1(#2); l2=data$2(#2); l3=data$3(#2);  u=data$4(#3);"
		  "t=comp(Base(#1).Grad(#1).Base(#2).Base(#3));"
		  "M$1(#1,#1)+=t(:,:,1,i,p).l1(i).u(p)+t(:,:,2,i,p).l2(i).u(p)+t(:,:,3,i,p).l3(i).u(p);"); 
	assem.push_mi(mim);
	assem.push_mf(mf_c);
	assem.push_mf(mf_data);
	assem.push_mf(mf_u);
	assem.push_data(lambdax);
	assem.push_data(lambday);
	assem.push_data(lambdaz);
	assem.push_data(U);
	assem.push_mat(D);
	assem.assembly(rg);

}


/*! Build the mixed boundary conditions (weak form) and dirichlet (strong form) for vessels
    @f$ M=\int_{\mathcal{E}_{MIX}} \beta~c~b~d\sigma@f$ and
    @f$ F=\int_{\mathcal{E}_{MIX}} \beta~c0~b~d\sigma@f$
 */
/*!
	@param M        BC contribution to mass matrix
	@param F        BC contribution to rhs
	@param mim      The integration method to be used
	@param mf_c     The finite element method for the concentration @f$c@f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of network boundary points
	@param beta     The beta value for mix condition @f$c0@f$

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void
asm_network_bc_transp
	(MAT & M, VEC & F,
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & BC,
	 const scalar_type beta) 
{ 
	GMM_ASSERT1(mf_c.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
	GMM_ASSERT1(mf_data.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");


	for (size_type bc=0; bc < BC.size(); bc++) { 
		GMM_ASSERT1(mf_c.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC[bc].label=="DIR") { // Dirichlet BC
			// Add cv_in contribution to Fv, and add a penalty as a mass matrix
			VEC BC_temp(mf_data.nb_dof(), BC[bc].value);
			getfem::assembling_Dirichlet_condition(M, F, mf_c, BC[bc].rg, BC_temp);
			gmm::clear(BC_temp);				
		} 
		else if (BC[bc].label=="MIX") { // Robin BC
			VEC BETA(mf_data.nb_dof(), beta);
			getfem::asm_mass_matrix_param(M, mim, mf_c, mf_data, BETA,mf_c.linked_mesh().region(BC[bc].rg) );
			
			VEC BETA_C0(mf_data.nb_dof(), beta*BC[bc].value);
			asm_source_term(F,mim, mf_c, mf_data,BETA_C0);
	
		}
		else if (BC[bc].label=="INT") { // Internal Node
			GMM_WARNING1("internal node passed as boundary.");
		}
		else if (BC[bc].label=="JUN") { // Junction Node
			GMM_WARNING1("junction node passed as boundary.");
		}
		else {
			GMM_ASSERT1(0, "Unknown Boundary Condition"<< BC[bc].label << endl);
		}
	}

}


} /* end of namespace */

#endif
