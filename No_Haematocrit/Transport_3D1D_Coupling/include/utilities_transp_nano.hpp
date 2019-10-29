/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
 "Computational models for nanoparticle transport in the vascular system"
    	       Master thesis in Mathematical Engineering
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
         Copyright (C) 2017 Annagiulia Tiozzo, Federica Laurino
======================================================================*/
/*! 
  @file   utilities_transp_nano.hpp
  @author Annagiulia Tiozzo <annagiulia92t@gmail.com>
	  Federica Laurino <federica.laurino@polimi.it>
  @date   April 2017.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */
#ifndef M3D1D_UTILITIES_TRANSP_NANO_HPP_
#define M3D1D_UTILITIES_TRANSP_NANO_HPP_

#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_mesh.h>
#include <gmm/gmm.h>
#include <math.h>

namespace getfem {
//! Compute the Wall Shear Stress: 	
//! @f$ W = 4 vel / Rad @f$
/*
 @param WSR  wall shear stress
 @param vel velocity 
 @param Rad radius of the branch
*/

template <typename VEC>
void wall_shear_stress	
	(VEC & W,
	const scalar_type & Rad,
	VEC & vel
	)
{	
	int lengthW = gmm::mat_nrows(gmm::col_vector(W));
	int lengthvel= gmm::mat_nrows(gmm::col_vector(vel));

	gmm::add(vel, W);
	gmm::scale(W, 4.0/Rad); 
};


//Compute the probability of adhesion P_a (s) from the exponential function	
//! @f$ P=ml~Ka~mr~\pi~r0^2~ \exp(-beta~mu~WSR/alpha) @f$
/*
 @param P 	The probability of adhesion
 @param WSR  	The wall shear rate
 @param ml 	The ligand density
 @param Ka 	The affinity constant of the ligand receptor bond
 @param alpha 	Coefficient for mr
 @param r0	The radius of the adhesion point
 @param beta	Coefficient for @f$ \frac{6F\lambda}{k_B T} @f$
 @param mu	The viscosity
 @param mr      The receptor density
*/

template <typename VEC>
void probability_adhesion0
	(VEC & P,
	VEC & W,
	const scalar_type & ml,
	const scalar_type & Ka,
	const scalar_type & alpha,
	const scalar_type & r0,
	const scalar_type & beta,
	const scalar_type & mu,
	const scalar_type & mr
	)

{ 	int lengthP= gmm::mat_nrows(gmm::col_vector(P));
	int lengthW= gmm::mat_nrows(gmm::col_vector(W));
	GMM_ASSERT1(lengthP==lengthW, "impossible to operate on the vectors: different lengths");

	vector_type temporary(lengthP); gmm::clear(temporary); 
	vector_type temporary2(lengthP); gmm::clear(temporary2); 
	gmm::add(W,temporary);

	gmm::scale(temporary, -(beta*mu/alpha));

	for(int i=0; i<lengthP; i++){
	temporary2[i]=exp(temporary[i]);
	}


	gmm::scale(temporary2, (ml*Ka*mr*pi*pow(r0,2.0)));
	gmm::add(temporary2,P);
	gmm::clear(temporary);
	gmm::clear(temporary2);	


};

//Compute the probability of adhesion P_a (s) from the Lattice Boltzmann approach 	
/*
 @param P 	probability of adhesion
 @param Re  	Reynolds number
 @param Pa_min 	P_a(Re_min)
 @param Pa_med 	P_a(Re_med)
 @param Pa_max 	P_a(Re_max)
*/

template <typename VEC>
void probability_adhesion
	(VEC & P,
	VEC & Re,
	scalar_type Pa_min,
	scalar_type Pa_med,
	scalar_type Pa_max
	)

{ 	int lengthP= gmm::mat_nrows(gmm::col_vector(P));
	int lengthRe= gmm::mat_nrows(gmm::col_vector(Re));
	GMM_ASSERT1(lengthP==lengthRe, "impossible to operate on the vectors: different lengths");

	scalar_type Re_min=0.01;
	scalar_type Re_med=0.1;
	scalar_type Re_max=1.0;


	for(int i=0; i<lengthP; i++){
		if(Re[i]<=Re_min)
			P[i]=1.0; 
		if(Re[i]<=Re_med && Re[i]>Re_min)
			P[i]=(Pa_med-Pa_min)/(Re_med-Re_min)*(Re[i]-Re_min)+Pa_min;
		if(Re[i]<=Re_max && Re[i]>Re_med)
			P[i]=(Pa_max-Pa_med)/(Re_max-Re_med)*(Re[i]-Re_med)+Pa_med;
		if(Re[i]>Re_max)
			P[i]=0.0; 
	}




};

//Compute the Reynolds number: 	
//! @f$ Re = vel \rho Rad / \mu @f$

/*
 @param Re 	Reynolds number
 @param vel  	Velocity
 @param mu 	Viscosity
 @param rho 	Density
*/
template <typename VEC>
void reynolds
	(VEC & Re,
	 VEC & vel,
	const scalar_type & Rad,
	const scalar_type & mu,
	const scalar_type & rho
	)
{
	int lengthRe= gmm::mat_nrows(gmm::col_vector(Re));
	int lengthVel= gmm::mat_nrows(gmm::col_vector(vel));

	gmm::add(vel, Re);
	gmm::scale(Re, rho*Rad/mu);


};


} /* end of namespace */



#endif
