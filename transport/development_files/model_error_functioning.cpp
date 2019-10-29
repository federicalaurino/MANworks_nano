/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2018 Stefano Brambilla
======================================================================*/
/*! 
  @file   model_error.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   July 2018
  @brief  Methods of the main class for computing model error.
 */
 
 #include <transport3d1d.hpp>
 #include <AMG_Interface.hpp>
 #include <cmath>
 #include "gmm/gmm_inoutput.h"
 #include "getfem/getfem_import.h"


 namespace getfem {


//////////////////////////////////////////////////
// Some global variables
const int A0=0;
const int A1=1;
const int A2=2;
const int A3=3;

const int PRIMAL=1;
const int DUAL=2;

///////////////////////////////////////////
// declaration of some useful parameters for exact solution, source terms, etc
double CC=1.0, kk=1.0, RR=1.0;
std::string expr="";

// Function declaration for source terms, etc

/* In this functions, we use the useful muParser library.
 * For further informations, check the following link: 
 * http://beltoforion.de/article.php?a=muparser
 * Anyway, you can pass whatever expression do you need, 
 * with a wide range of possiblities, such as:
 * sin,cos,tan,log,ln,log2,log10,exp,sqrt,sign,abs,
 * rint (round to integer),min,max,sum,avg,
 * =,==,&&,||,<=,>=,+,-,*,/,^.
*/

//! Function for source term in vessels
double g_function(const bgeot::base_node & x){

	//Define p
	mu::Parser p;
	//Define variables
		/* We define both coordinates x,y,z,
		 * which are needed for the function,
		 * and other useful parameters, namely C,R,k,
		 * that you can pass in the .param file.
		 * Note that you can add in this file any other 
		 * variable you want. */ 
	double xx=x[0], yy=x[1], zz=x[2];
	p.DefineVar("x", &xx);
	p.DefineVar("y", &yy);
	p.DefineVar("z", &zz);
	p.DefineVar("C", &CC);	
	p.DefineVar("k", &kk);
	p.DefineVar("R", &RR);
	//p.DefineVar("new variable", double variable);	
	
	//Define the expression
		/* The string should be already modified
		 * in the code */
	p.SetExpr(expr);
	return p.Eval();
}

// Used in model_error:
//! Function for source term in tissue 
double f_function(const bgeot::base_node & x){

	//Define p
	mu::Parser p;
	//Define variables
		/* We define both coordinates x,y,z,
		 * which are needed for the function,
		 * and other useful parameters, namely C,R,k,
		 * that you can pass in the .param file.
		 * Note that you can add in this file any other 
		 * variable you want. */ 
	double xx=x[0], yy=x[1], zz=x[2];
	p.DefineVar("x", &xx);
	p.DefineVar("y", &yy);
	p.DefineVar("z", &zz);
	p.DefineVar("C", &CC);	
	p.DefineVar("k", &kk);
	p.DefineVar("R", &RR);
	//p.DefineVar("new variable", double variable);	
	
	//Define the expression
		/* The string should be already modified
		 * in the code */
	p.SetExpr(expr);
	return p.Eval();
}

	/////////////////////////////////////////////////////
	// Methods for computing reduced (primal and dual) solution

	//! Initialize problem U1
	void transport3d1d::init_model(int argc, char *argv[]){

	PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	import_data_transp();
	//2. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh_transp();
	//3. Set finite elements and integration methods
	set_im_and_fem_model();
 	//4. Build problem parameters
	build_param_transp();
	//5. Build the list of tissue boundary data
	build_tissue_boundary_transp();
	//6. Build the list of tissue boundary (and junction) data
 	build_vessel_boundary_transp();
};
	//3. Set finite elements and integration methods
	void transport3d1d::set_im_and_fem_model(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	
	
	pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
	pfem pf_Cv = fem_descriptor(descr_transp.FEM_TYPEV_C); 

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for tissue ..." << endl;
	#endif
		

	//Define mf_Ct_Omega and mf_Ct_Sigma
	mf_Ct.set_finite_element(mesht.convex_index(), pf_Ct);
	mf_Ct_Omega.set_finite_element(mesht.region(descr_transp.OMEGA).index(), pf_Ct);
	mf_Ct_Sigma.set_finite_element(mesht.region(descr_transp.SIGMA).index(), pf_Ct); 

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif

	mf_Cv.set_finite_element(meshv.convex_index(), pf_Cv);

	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif

	dof_transp.set(mf_Ct, mf_Cv, mf_Ct_Omega, mf_Ct_Sigma);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof_transp;
	#endif
	
	mimv.clear();
	mf_Uvi.clear();
	mf_Pv.clear();
	mf_coefv.clear();
	mf_coefvi.clear();

	mimt.clear();
	mf_Ut.clear();
	mf_Pt.clear();
	mf_coeft.clear();
	problem3d1d::set_im_and_fem();


};
	//! Assemble problem U1
	void transport3d1d::assembly_model(const size_type ASSUMPTION){

	// define dofs
	size_type dof_t;
	size_type dof_v;

	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		dof_t  = dof_transp.Ct_Omega();	}
	else{
		dof_t  = dof_transp.Ct();	}

 	if(ASSUMPTION==A0){
		dof_v  = dof_transp.Ct_Sigma();	}
	else{
		dof_v  = dof_transp.Cv();	}

	size_type dof_tot= dof_t+dof_v;

	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM_transp, UM_transp, FM_transp ..." << endl;
	#endif
	gmm::resize(AM_transp, dof_tot, dof_tot);	gmm::clear(AM_transp);
	gmm::resize(UM_transp, dof_tot);		gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_tot);	 	gmm::clear(FM_transp);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM_transp ..." << endl;
	#endif
	// Mass(time derivative)  matrix for the interstitial problem
	sparse_matrix_type At(dof_t, dof_t);gmm::clear(At);
	// Mass (time derivative)  matrix for the network problem
	sparse_matrix_type Av(dof_v, dof_v);gmm::clear(Av);	

	// Tissue-to-tissue exchange matrix on Gamma
	sparse_matrix_type Mtt(dof_t, dof_t);gmm::clear(Mtt);
	// Tissue-to-vessel exchange matrix on Gamma
	sparse_matrix_type Mtv(dof_t, dof_v);gmm::clear(Mtv);
	// Vessel-to-tissue exchange matrix on Gamma
	sparse_matrix_type Mvt(dof_v, dof_t);gmm::clear(Mtv);
	// Vessel-to-vessel exchange matrix on Gamma
	sparse_matrix_type Mvv(dof_v, dof_v);gmm::clear(Mvv);

	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof_t, dof_t);gmm::clear(Btt);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof_t, dof_v);gmm::clear(Btv);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_v, dof_t);gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof_v, dof_v);gmm::clear(Bvv);

	// Aux tissue-to-vessel averaging matrix (Circumference)
	sparse_matrix_type Mbar(dof_v, dof_t);gmm::clear(Mbar);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof_v, dof_t);gmm::clear(Mlin);
	// Aux tissue-to-vessel averaging matrix (Section)
	sparse_matrix_type Mbarbar(dof_v, dof_transp.Ct());gmm::clear(Mbarbar); //Whatever the case, this matrix should be defined on the whole 3d mesh
	// Aux tissue source vector
	vector_type Ft(dof_t); gmm::clear(Ft);
	// Aux vessels source vector
	vector_type Fv(dof_v); gmm::clear(Fv);



/*
	// Tissue-to-tissue exchange matrix on Gamma
	sparse_matrix_type M_gamma_singolo(dof_transp.Ct(), dof_transp.Ct());gmm::clear(M_gamma_singolo);
	sparse_matrix_type M_gamma_doppio(dof_transp.Ct(), dof_transp.Ct());gmm::clear(M_gamma_doppio);
	sparse_matrix_type M_gamma_vuoto(dof_transp.Ct_Omega(), dof_transp.Ct_Omega());gmm::clear(M_gamma_vuoto);
	sparse_matrix_type M_gamma_vuoto_doppio(dof_transp.Ct_Omega(), dof_transp.Ct_Omega());gmm::clear(M_gamma_vuoto_doppio);
	getfem::asm_mass_matrix (M_gamma_singolo, mimt, mf_Ct, descr_transp.GAMMA);
	getfem::asm_mass_matrix (M_gamma_doppio, mimt, mf_Ct, descr_transp.GAMMA+1);
	getfem::asm_mass_matrix (M_gamma_vuoto, mimt, mf_Ct_Omega, descr_transp.GAMMA);
	getfem::asm_mass_matrix (M_gamma_vuoto_doppio, mimt, mf_Ct_Omega, descr_transp.GAMMA+1);
	gmm::MatrixMarket_IO::write("M_gamma_singolo.mm" , M_gamma_singolo);
	gmm::MatrixMarket_IO::write("M_gamma_doppio.mm" , M_gamma_doppio);
	gmm::MatrixMarket_IO::write("M_gamma_vuoto.mm" , M_gamma_vuoto);
	gmm::MatrixMarket_IO::write("M_gamma_vuoto_doppio.mm" , M_gamma_vuoto_doppio);
	vector_type Grad_gamma_singolo(dof_transp.Ct());gmm::clear(Grad_gamma_singolo);
	vector_type Grad_gamma_doppio(dof_transp.Ct());gmm::clear(Grad_gamma_doppio);
	vector_type Grad_gamma_vuoto(dof_transp.Ct_Omega());gmm::clear(Grad_gamma_vuoto);
	vector_type Grad_gamma_vuoto_doppio(dof_transp.Ct_Omega());gmm::clear(Grad_gamma_vuoto_doppio);
	vector_type Grad_faccia0(dof_transp.Ct());gmm::clear(Grad_faccia0);
	vector_type ones(3,1);
	getfem::asm_homogeneous_normal_source_term (Grad_gamma_singolo, mimt, mf_Ct, ones, descr_transp.GAMMA);
	getfem::asm_homogeneous_normal_source_term (Grad_gamma_doppio, mimt, mf_Ct, ones, descr_transp.GAMMA+1);
	getfem::asm_homogeneous_normal_source_term (Grad_gamma_vuoto, mimt, mf_Ct_Omega, ones, descr_transp.GAMMA);
	getfem::asm_homogeneous_normal_source_term (Grad_gamma_vuoto_doppio, mimt, mf_Ct_Omega, ones, descr_transp.GAMMA+1);
	getfem::asm_homogeneous_normal_source_term (Grad_faccia0, mimt, mf_Ct, ones,0);
	std::ofstream out1(descr_transp.OUTPUT+"Grad_gamma_singolo.txt");
		out1 << gmm::col_vector(Grad_gamma_singolo);
		out1.close(); 
	std::ofstream out2(descr_transp.OUTPUT+"Grad_gamma_doppio.txt");
		out2 << gmm::col_vector(Grad_gamma_doppio);
		out2.close(); 
	std::ofstream out3(descr_transp.OUTPUT+"Grad_gamma_vuoto.txt");
		out3 << gmm::col_vector(Grad_gamma_vuoto);
		out3.close(); 
	std::ofstream out4(descr_transp.OUTPUT+"Grad_gamma_vuoto_doppio.txt");
		out4 << gmm::col_vector(Grad_gamma_vuoto_doppio);
		out4.close(); 
	std::ofstream out5(descr_transp.OUTPUT+"Grad_faccia0.txt");
		out5 << gmm::col_vector(Grad_faccia0);
		out5.close(); 
*/


	//Assemble Perimeter, Area, and build Kt and Kv

	//Permeability
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	vector_type Kv(dof.coefv(),k);
	vector_type Kt(dof.coeft(),k);

	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);

	//Area = pi*R^2
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); 

 
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix At ..." << endl;
	#endif	
	// Assemble At, stiffness matrix for laplacian in tissue
	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct_Omega); 	}
	else{
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct      );    }
	
	gmm::add(At,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 


	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix Av ..." << endl;
	#endif	
	// Assemble Av, stiffness matrix for laplacian in vessel
	if(ASSUMPTION==A0){
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(Av,mimt,mf_Ct_Sigma); 	}
	else{
		getfem::asm_stiffness_matrix_for_laplacian(Av,mimv,mf_Cv,mf_coefv,Area);     }

	
	gmm::add(Av,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_t, dof_v), 
					gmm::sub_interval(dof_t, dof_v))); 


	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	if(ASSUMPTION==A0){

		//The mass terms on Gamma are easy to compute
		getfem::asm_mass_matrix_param (Btt, mimt, mf_Ct_Omega,mf_Ct_Omega, mf_coeft, Kt, descr_transp.GAMMA);
		getfem::asm_mass_matrix_param (Bvv, mimt, mf_Ct_Sigma,mf_Ct_Sigma, mf_coeft, Kt, descr_transp.GAMMA+1);

		//getfem::asm_mass_matrix_param (Btv, mimt, mf_Ct_Omega,mf_Ct_Sigma, mf_coeft, Kt );
		//getfem::asm_mass_matrix_param (Bvt, mimt, mf_Ct_Sigma,mf_Ct_Omega, mf_coeft, Kt);
	
	/* For exchange terms on Gamma, I must use interpolation matrix between Omega and Sigma.
	 * In fact, in the mesht, are stored only the elements, not their faces.
	 * i.e. elements from 1 to 10.000, and faces from 0 to 3 of each element. (there is not a unic ID for a face)
	 * Therefore, when i build Gamma to be the outer faces of Sigma, it takes the faces of the elements stored in Sigma. 
	 * If i buildt Gamma to be the outer faces of Omega (beware, you should take care of the exterior faces of the cube),
	 * this region would contain exactly the same faces, but with different ID.
	 * i.e. face 2 of element 120 in Omega == face 1 of element 6300 in Sigma, but GetFEM DOESN'T KNOW THAT!
	 * SOLUTION: I build region Gamma+1 to be the outer faces of Sigma AND outer faces of Omega.
	 * I build the mesh_fem of Gamma+1, and build the mass matrix on Gamma (M_Gamma).
	 * I assemble the interpolation matrixes from Gamma to Sigma and from Gamma to Omega (INT_Gamma_Sigma and INT_Gamma_Omega).
	 * (u_omega, v_sigma)_Gamma = INT_Gamma_Omega* (INT_Gamma_Sigma * M_Gamma)^T = Btv = (Bvt)^T 
	 */
		pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
		mesh_fem mf_Ct_Gamma(mesht);
		mf_Ct_Gamma.set_finite_element(mesht.region(descr_transp.GAMMA+1).index(), pf_Ct);
		/*mesh_fem mf_Ct_Gamma0(mesht);
		mf_Ct_Gamma0.set_finite_element(mesht.region(descr_transp.GAMMA).index(), pf_Ct);
cout <<"Gamma ha "<<mf_Ct_Gamma0.nb_dof()<<" dof e Gamma_doppio ha "<<mf_Ct_Gamma.nb_dof()<<" dof."<<endl;
cout <<"Omega ha "<<mf_Ct_Omega.nb_dof()<<" dof e Sigma ha "<<mf_Ct_Sigma.nb_dof()<<" dof."<<endl;
cout <<"Ct ha "<<mf_Ct.nb_dof()<<" dof e Cv ha "<<mf_Cv.nb_dof()<<" dof."<<endl;
*/
		sparse_matrix_type M_Gamma(mf_Ct_Gamma.nb_dof(),mf_Ct_Gamma.nb_dof()); gmm::clear(M_Gamma);
		sparse_matrix_type INT_Gamma_Sigma(mf_Ct_Sigma.nb_dof(),mf_Ct_Gamma.nb_dof()); gmm::clear(INT_Gamma_Sigma);
		sparse_matrix_type INT_Gamma_Omega(mf_Ct_Omega.nb_dof(),mf_Ct_Gamma.nb_dof()); gmm::clear(INT_Gamma_Omega);
		getfem::interpolation(mf_Ct_Gamma, mf_Ct_Sigma, INT_Gamma_Sigma);
		getfem::interpolation(mf_Ct_Gamma, mf_Ct_Omega, INT_Gamma_Omega);
			
		getfem::asm_mass_matrix(M_Gamma, mimt, mf_Ct_Gamma);	
		gmm::mult(INT_Gamma_Sigma, M_Gamma, INT_Gamma_Sigma); // M1 * M2 ---> M3			
		gmm::mult(INT_Gamma_Omega, gmm::transposed(INT_Gamma_Sigma), Btv);
		gmm::copy(gmm::transposed(Btv), Bvt); 

	}
	else if(ASSUMPTION==A1){

		asm_exchange_aux_mat_transp(Mbar, Mlin, 
			mimv, mf_Ct_Omega, mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);

	// What should i do if I want to save this Mbar (on Omega) and use it afterwards on whole mf_Ct?? 
	
/*	cout <<"Mbar e Mlin sono salvati, le dimensioni di Mbar sono: colonne = "<<gmm::mat_ncols(Mbar)<<" e righe = " <<gmm::mat_nrows(Mbar)<<endl; 
	cout <<"i gradi di liberta sono: mf_Cv = "<<dof_v<<" e mf_Ct = "<<dof_t<<endl;
	cout <<"max_norm(Mbar) = "<<gmm::mat_maxnorm(Mbar)<<"max_norm(Mlin) = "<<gmm::mat_maxnorm(Mlin)<<endl;
	cout <<"frob_norm(Mbar) = "<<gmm::mat_euclidean_norm(Mbar)<<"frob_norm(Mlin) = "<<gmm::mat_euclidean_norm(Mlin)<<endl;
	cout <<"mat_norm1(Mbar) = "<<gmm::mat_norm1(Mbar)<<"mat_norm1(Mlin) = "<<gmm::mat_norm1(Mlin)<<endl;	
	cout <<"mat_norminf(Mbar) = "<<gmm::mat_norminf(Mbar)<<"mat_norminf(Mlin) = "<<gmm::mat_norminf(Mlin)<<endl;
	
	sparse_matrix_type Momega(dof_transp.Ct(),dof_transp.Ct_Omega()); gmm::clear(Momega);
	getfem::interpolation(mf_Ct_Omega, mf_Ct, Momega, 2);
	gmm::resize(MLIN, dof_transp.Cv(), dof_transp.Ct());
	gmm::resize(MBAR, dof_transp.Cv(), dof_transp.Ct());
	gmm::mult(Momega, Mlin, MLIN);
	gmm::mult(Momega, Mbar, MBAR); 
	cout <<endl<<"MBAR e MLIN sono salvati, le dimensioni di MBAR sono: colonne = "<<gmm::mat_ncols(MBAR)<<" e righe = " <<gmm::mat_nrows(MBAR)<<endl; 
	cout <<"i gradi di liberta sono: mf_Cv = "<<dof_transp.Cv()<<" e mf_Ct = "<<dof_transp.Ct()<<endl;
	cout <<"max_norm(MBAR) = "<<gmm::mat_maxnorm(MBAR)<<"max_norm(MLIN) = "<<gmm::mat_maxnorm(MLIN)<<endl;
	cout <<"frob_norm(MBAR) = "<<gmm::mat_euclidean_norm(MBAR)<<"frob_norm(MLIN) = "<<gmm::mat_euclidean_norm(MLIN)<<endl;
	cout <<"mat_norm1(MBAR) = "<<gmm::mat_norm1(MBAR)<<"mat_norm1(MLIN) = "<<gmm::mat_norm1(MLIN)<<endl;	
	cout <<"mat_norminf(MBAR) = "<<gmm::mat_norminf(MBAR)<<"mat_norminf(MLIN) = "<<gmm::mat_norminf(MLIN)<<endl;
*/

/*
gmm::mat_euclidean_norm(M) // Euclidean norm of matrix ‘‘M‘‘
// (called also Frobenius norm).
gmm::mat_maxnorm(M) // Max norm (defined as max(|m_ij|; i,j ))
gmm::mat_norm1(M)
// max(sum(|m_ij|, i), j)
gmm::mat_norminf(M) // max(sum(|m_ij|, j), i)
*/

/*
	// Aux tissue-to-vessel averaging matrix (Circumference)
	sparse_matrix_type Mbar2(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mbar2);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin2(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mlin2);
	asm_exchange_aux_mat_transp(Mbar2, Mlin2, 
			mimv, mf_Ct,       mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);
		
	cout <<endl<<"MBAR2 e MLIN2 sono salvati, le dimensioni di MBAR2 sono: colonne = "<<gmm::mat_ncols(Mbar)<<" e righe = " <<gmm::mat_nrows(Mbar)<<endl; 
	cout <<"i gradi di liberta sono: mf_Cv = "<<dof_transp.Cv()<<" e mf_Ct = "<<dof_transp.Ct()<<endl;
	cout <<"max_norm(MBAR2) = "<<gmm::mat_maxnorm(Mbar2)<<"max_norm(MLIN2) = "<<gmm::mat_maxnorm(Mlin2)<<endl<<endl;
	cout <<"frob_norm(MBAR2) = "<<gmm::mat_euclidean_norm(Mbar2)<<"frob_norm(MLIN2) = "<<gmm::mat_euclidean_norm(Mlin2)<<endl;
	cout <<"mat_norm1(MBAR2) = "<<gmm::mat_norm1(Mbar2)<<"mat_norm1(MLIN2) = "<<gmm::mat_norm1(Mlin2)<<endl;	
	cout <<"mat_norminf(MBAR2) = "<<gmm::mat_norminf(Mbar2)<<"mat_norminf(MLIN2) = "<<gmm::mat_norminf(Mlin2)<<endl;
	cout<<"Calcolo MBAR-MBAR2"<<endl;
	
	gmm::add(gmm::scaled(MLIN,-1), Mlin2); // V1 + V2 --> V2	
	gmm::add(gmm::scaled(MBAR,-1), Mbar2); 
	cout <<"max_norm(MBAR - MBAR2) = "<<gmm::mat_maxnorm(Mbar2)<<"max_norm(MLIN - MLIN2) = "<<gmm::mat_maxnorm(Mlin2)<<endl<<endl;
	cout <<"frob_norm(MBAR - MBAR2) = "<<gmm::mat_euclidean_norm(Mbar2)<<"frob_norm(MLIN - MLIN2) = "<<gmm::mat_euclidean_norm(Mlin2)<<endl;
	cout <<"mat_norm1(MBAR - MBAR2) = "<<gmm::mat_norm1(Mbar2)<<"mat_norm1(MLIN - MLIN2) = "<<gmm::mat_norm1(Mlin2)<<endl;	
	cout <<"mat_norminf(MBAR - MBAR2) = "<<gmm::mat_norminf(Mbar2)<<"mat_norminf(MLIN - MLIN2) = "<<gmm::mat_norminf(Mlin2)<<endl;
*/
	}
	else{
		if( (gmm::mat_nrows(MBAR)==dof_v) && (gmm::mat_ncols(MBAR)==dof_t) 
		&&  (gmm::mat_nrows(MLIN)==dof_v) && (gmm::mat_ncols(MLIN)==dof_t)
		&&  (gmm::mat_maxnorm(MBAR)>=1e-16) && (gmm::mat_maxnorm(MLIN)>=1e-16) ) //Mbar and Mlin are already defined
			{gmm::copy(MBAR, Mbar); gmm::copy(MLIN, Mlin); }
		else{ //Mbar and Mlin are not already defined
			asm_exchange_aux_mat_transp(Mbar, Mlin, 
				mimv, mf_Ct,       mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);
			gmm::resize(MLIN, dof_transp.Cv(), dof_transp.Ct());
			gmm::resize(MBAR, dof_transp.Cv(), dof_transp.Ct());
			gmm::copy(Mbar,MBAR); gmm::copy(Mlin,MLIN);
		}
	}


	bool NEWFORM = 1;//PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");	
	if(ASSUMPTION != A0){									   
		asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, Mbar, Mlin, Perimeter, NEWFORM);}
	// Copying Bvt
	gmm::add(scaled(Bvt,-1),								
			  gmm::sub_matrix(AM_transp, 
			  		gmm::sub_interval(dof_t, dof_v),
					gmm::sub_interval(0, dof_t)));
	// Copying Bvv
	gmm::add(Bvv,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_t, dof_v), 
					gmm::sub_interval(dof_t, dof_v))); 
	if(ASSUMPTION==A1){
		gmm::resize(BVT_Omega, dof_transp.Cv(), dof_transp.Ct_Omega());
		gmm::copy(Bvt, BVT_Omega);}
/*
	// Print the first elements of Bvv, Bvt, Btv and Btt (check if they are all-zero matrixes)
	int cont_max=0;
	cout<<"Ecco i primi 5 elementi non nulli di Btt: "<<endl;
	for (int i1=0; (i1<dof_t)&&(cont_max<5); i1++){
		for (int i2=0; (i2<dof_t)&&(cont_max<5); i2++){
	if(Btt(i1,i2)!=0){ cout<<"elemento ( "<<i1<<" , "<<i2<<" ) è uguale a: "<<Btt(i1,i2)<<endl; cont_max++;}
	}}
	cont_max=0;
	cout<<"Ecco i primi 5 elementi non nulli di Btv: "<<endl;
	for (int i1=0; (i1<dof_t)&&(cont_max<5); i1++){
		for (int i2=0; (i2<dof_v)&&(cont_max<5); i2++){
	if(Btv(i1,i2)!=0){ cout<<"elemento ( "<<i1<<" , "<<i2<<" ) è uguale a: "<<Btv(i1,i2)<<endl; cont_max++;}
	}}
	cont_max=0;
	cout<<"Ecco i primi 5 elementi non nulli di Bvt: "<<endl;
	for (int i1=0; (i1<dof_v)&&(cont_max<5); i1++){
		for (int i2=0; (i2<dof_t)&&(cont_max<5); i2++){
	if(Bvt(i1,i2)!=0){ cout<<"elemento ( "<<i1<<" , "<<i2<<" ) è uguale a: "<<Bvt(i1,i2)<<endl; cont_max++;}
	}}
	cont_max=0;
	cout<<"Ecco i primi 5 elementi non nulli di Bvv: "<<endl;
	for (int i1=0; (i1<dof_v)&&(cont_max<5); i1++){
		for (int i2=0; (i2<dof_v)&&(cont_max<5); i2++){
	if(Bvv(i1,i2)!=0){ cout<<"elemento ( "<<i1<<" , "<<i2<<" ) è uguale a: "<<Bvv(i1,i2)<<endl; cont_max++;}
	}}
*/
	if(ASSUMPTION==A1 || ASSUMPTION==A2){
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling mass matrices on Gamma..." << endl;
	#endif
		// Assemble Mtt, mass matrix on tissue
		if(ASSUMPTION==A1){
			getfem::asm_mass_matrix_param (Mtt, mimt, mf_Ct_Omega, mf_coeft, Kt, descr_transp.GAMMA);}
		else{
			getfem::asm_mass_matrix_param (Mtt, mimt, mf_Ct,       mf_coeft, Kt, descr_transp.GAMMA);}

		// Assemble Mtv, exchange matrix on gamma
		getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Perimeter);
		gmm::mult(gmm::transposed(Mbar), Mvv, Mtv); // M1 * M2 ---> M3	
	
		// Copying Mtt
		gmm::add(Mtt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 
		// Copying Mtv

		gmm::add(scaled(Mtv,-1),
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t),
					gmm::sub_interval(dof_t, dof_v))); 

	} else { //ASSUMPTION==A3 || ASSUMPTION==A0
	// Copying Btt
		gmm::add(Btt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 

		if(ASSUMPTION==A3){
			gmm::resize(BTT, dof_transp.Ct(), dof_transp.Ct());	gmm::clear(BTT);
			gmm::copy(Btt,BTT);
		}

		// Copying Btv

		gmm::add(scaled(Btv,-1),
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t),
					gmm::sub_interval(dof_t, dof_v))); 

	}

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling source terms ..." << endl;
	#endif
	// Update potential paramters in F and G
	kk=k;
	CC = PARAM.real_value("Cv", "value of concentration in the vessel");
	RR = param.R(0);

	//// Assemble F: source term in tissue
	vector_type F(dof.coeft());
	//Vector F is buildt with f_function via muparser
	expr = PARAM.string_value("F_EXPR", "Expression for source term in tissue (muparser)");
	interpolation_function(mf_coeft, F, f_function ); 
	//Vector F is assembled on Ct_Omega or on whole Ct?	
	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		getfem::asm_source_term(Ft, mimt, mf_Ct_Omega, mf_coeft, F);}
	else{
		getfem::asm_source_term(Ft, mimt, mf_Ct,       mf_coeft, F);}	
	

	//// Assemble G: source term in vessel
	if(ASSUMPTION==A0){	// G is defined on Sigma
		vector_type G(dof.coeft());
		//Vector G is buildt with g_function via muparser
		expr = PARAM.string_value("G_EXPR", "Expression for source term in vessel (muparser)");
		interpolation_function(mf_coeft, G, g_function );  

		//Assemble source term on Sigma
		getfem::asm_source_term(Fv, mimt, mf_Ct_Sigma, mf_coeft, G);

		}
	else{			// G is defined on Lambda

		vector_type G(dof_transp.Cv());gmm::clear(G);
		expr = PARAM.string_value("G_EXPR", "Expression for source term in vessel (muparser)");
		bool G_CONSTANT = PARAM.int_value("G_CONSTANT", "flag for using constant source term g in vessel");

		// If G is constant in the cross section, i just interpolate the g_function in Lambda
		if(G_CONSTANT){
			interpolation_function(mf_Cv, G, g_function );	
		}
		else{	//If G is not constant, I must build Mbarbar
			vector_type Gt(dof_transp.Ct()); gmm::clear(Gt);	
			interpolation_function(mf_Ct, Gt, g_function ); 		
			
			// Build Mbarbar 
			if(gmm::mat_nrows(MBARBAR)==dof_transp.Cv()
			&& gmm::mat_ncols(MBARBAR)==dof_transp.Ct() 
			&& gmm::mat_maxnorm(MBARBAR)>=1e-16) //Mbarbar is already defined
				{gmm::copy(MBARBAR, Mbarbar); }
			else { 				     //Mbarbar is not already defined

				gmm::resize(MBARBAR, dof_transp.Cv(), dof_transp.Ct());	
				gmm::clear(MBARBAR);
				asm_exchange_aux_mat_bar_bar(Mbarbar,mimt, mf_Ct, mf_Cv,mf_coefv, param.R(), descr.NInt, descr_transp.NIntA, nb_branches);
				gmm::copy(Mbarbar, MBARBAR);  
				}
			gmm::mult(Mbarbar, Gt, G);		
		}
		
		//Build source term
		vector_type Area_Cv(dof_transp.Cv());gmm::clear(Area_Cv);
		getfem::interpolation(mf_coefv, mf_Cv, Area, Area_Cv);
		gmm::vscale(Area_Cv, G); //G=G*pi*R^2
		getfem::asm_source_term(Fv, mimv, mf_Cv, mf_Cv, G);

	} 


	// Add source term to RHS
	gmm::add(Ft, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(0,dof_t)));

	gmm::add(Fv, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(dof_t, dof_v)));


	// De-allocate memory
	gmm::clear(At);	   gmm::clear(Av); 
	gmm::clear(Mtt);   gmm::clear(Mtv); 
	gmm::clear(Mvt);   gmm::clear(Mvv); 	    
	gmm::clear(Mbar);  gmm::clear(Mlin);   gmm::clear(Mbarbar);  
	gmm::clear(Btt);   gmm::clear(Btv);
	gmm::clear(Bvt);   gmm::clear(Bvv);

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling boundary conditions ..." << endl;
	#endif	
	gmm::resize(AM_temp, dof_tot,
			     dof_tot);	gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_tot); 	gmm::clear(FM_temp);
	gmm::copy(AM_transp,AM_temp);
	gmm::copy(FM_transp,FM_temp);

	//Boundary conditions: 
	//On vessel, we have homogeneous neumann conditions: we simply don't implement the boundary term due to diffusion.

	//Homogeneous Dirichlet on all the faces of tissue	
	vector_type Dirichlet_null(dof_t,0);
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {

		if(ASSUMPTION==A1 || ASSUMPTION==A0){
		getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct_Omega, BCt_transp[bc].rg, Dirichlet_null);	}
		else{
		getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct,       BCt_transp[bc].rg, Dirichlet_null);	}
	} 
	

};

	//! Assemble problem Z1
	void transport3d1d::assembly_dual_model(const size_type ASSUMPTION){

	// define dofs
	size_type dof_t;
	size_type dof_v;

	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		dof_t  = dof_transp.Ct_Omega();	}
	else{
		dof_t  = dof_transp.Ct();	}

 	if(ASSUMPTION==A0){
		dof_v  = dof_transp.Ct_Sigma();	}
	else{
		dof_v  = dof_transp.Cv();	}

	size_type dof_tot= dof_t+dof_v;

	gmm::clear(UM_transp);

	// Aux tissue source vector
	vector_type Jt(dof_t); gmm::clear(Jt);
	// Aux vessels source vector
	vector_type Jv(dof_v); gmm::clear(Jv);


	// prepare AM and Fm temp
	gmm::resize(AM_temp, dof_tot,
			     dof_tot);	gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_tot); 	gmm::clear(FM_temp);
	gmm::copy(AM_transp,AM_temp);
	gmm::copy(FM_transp,FM_temp);

	gmm::scale(gmm::sub_matrix(AM_temp, 
					gmm::sub_interval(0, dof_t),
					gmm::sub_interval(dof_t, dof_v))
		 ,0);
	gmm::scale(gmm::sub_matrix(AM_temp, 
			  		gmm::sub_interval(dof_t, dof_v),
					gmm::sub_interval(0, dof_t))
		 ,0);

		#ifdef M3D1D_VERBOSE_
	cout << "  Assembling source terms ..." << endl;
	#endif
	string functional= PARAM.string_value("FUNCTIONAL", "descriptor for the functional J of the error you wanto to estimate");

	if(functional=="MEAN_VALUE"){
	// Aux tissue source vector
	vector_type Ft(dof.coeft(),1);
	// Aux vessels source vector
	vector_type Fv(dof.coefv(),1);

	// Assemble Jt: source term in tissue
	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		getfem::asm_source_term(Jt,mimt, mf_Ct_Omega, mf_coeft, Ft);}
	else{
		getfem::asm_source_term(Jt,mimt, mf_Ct,       mf_coeft, Ft);}

	// Assemble Jv: source term in vessel
	if(ASSUMPTION==A0){
		getfem::asm_source_term(Jv,mimt, mf_Ct_Sigma, mf_coeft, Ft);}
	else{
		getfem::asm_source_term(Jv,mimv, mf_Cv, mf_coefv, Fv ); }

	}
	else if (functional=="MEAN_STRESS"){

	vector_type ones(3,1);
	// Assemble Jt: source term in tissue
	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		getfem::asm_homogeneous_normal_source_term (Jt, mimt, mf_Ct_Omega, ones, descr_transp.GAMMA);}
	else{
		getfem::asm_homogeneous_normal_source_term (Jt, mimt, mf_Ct      , ones, descr_transp.GAMMA);}

	// Assemble Jv: source term in vessel
	if(ASSUMPTION==A0){
		getfem::asm_homogeneous_normal_source_term (Jv, mimt, mf_Ct_Sigma, ones, descr_transp.GAMMA);}
	else{
		for (size_type bc=0; bc < BCv_transp.size(); bc++) { 
			getfem::asm_homogeneous_normal_source_term (Jv, mimv, mf_Cv, ones, BCv_transp[bc].rg);}} 

	};

	//Add to FM_temp

	gmm::add(Jt, 
			gmm::sub_vector(FM_temp,
					gmm::sub_interval(0,dof_t)));

	gmm::add(Jv, 
			gmm::sub_vector(FM_temp,
					gmm::sub_interval(dof_t, dof_v)));


	//Homogeneous Dirichlet on all the faces of tissue
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling boundary conditions ..." << endl;
	#endif	

	vector_type Dirichlet_null(dof_t,0);
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
	if(ASSUMPTION==A1){
	getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct_Omega, BCt_transp[bc].rg, Dirichlet_null);}
	else{
	getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct,       BCt_transp[bc].rg, Dirichlet_null);}
			
	} 

/*
	if(functional=="L2NORM")
	// Assemble F: source term in tissue
	sparse_matrix_type JJt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(JJt);
		
	generic_assembly L2normt;
	L2normt.push_mi(mimt);
	L2normt.push_mf(mf_Ct);
	L2normt.push_mat(JJt);
	L2normt.set("V(#1,#1)+=comp(Base(#1).Base(#1))(:,:)");
	L2normt.assembly();
	
	//gmm::add(gmm::mat_trace(JJt),Jt);
	
	for (int i=0; i>dof_transp.Ct(); i++){
		Jt[i]=JJt(i,i);
		Jt[i]=sqrt(Jt[i]);
	};
	// Assemble G: source term in vessel
	//! /todo Extend to variable source term: use muParser and build cross section average (operator Mbarbar)
	sparse_matrix_type JJv(dof_v, dof_v);gmm::clear(JJv);
		
	generic_assembly L2normv;
	L2normv.push_mi(mimv);
	L2normv.push_mf(mf_Cv);
	L2normv.push_mat(JJv);
	L2normv.set("V(#1,#1)+=comp(Base(#1).Base(#1))(:,:)");
	L2normv.assembly();
	
	//gmm::add(gmm::mat_trace(JJv),Jv);
	for (int i=0; i>dof_v; i++){
		Jv[i]=JJv(i,i);
		Jv[i]=sqrt(Jv[i]);
		
	};
*/	


};


	//! Solve problem U1
	bool transport3d1d::solve_model(const size_type ASSUMPTION,const  size_type VERSION){

  	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	double time = gmm::uclock_sec();
	gmm::clear(UM_transp);
	gmm::clean(AM_temp, 1E-12);
	gmm::clean(FM_temp, 1E-12);


	// define dofs
	size_type dof_t;
	size_type dof_v;

	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		dof_t  = dof_transp.Ct_Omega();	}
	else{
		dof_t  = dof_transp.Ct();	}
	
 	if(ASSUMPTION==A0){
		dof_v  = dof_transp.Ct_Sigma();	}
	else{
		dof_v  = dof_transp.Cv();	}

	size_type dof_tot= dof_t+dof_v;

	// Solve the system on AM_temp, UM_transp, FM_temp
	bool solved = solver(dof_t,dof_v);
	if (!solved) return false;
	//export solution
	#ifdef M3D1D_VERBOSE_
	std::cout<<"solved! going to export..."<<std::endl;
	#endif		



	switch(VERSION){

		case PRIMAL:{
		switch(ASSUMPTION){

			case A0: 
			{	gmm::resize(U0, dof_tot);	gmm::clear(U0);
				gmm::copy(UM_transp, U0);
				cout << endl<<"... time to solve U0: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A1: 
			{	gmm::resize(U1, dof_tot);	gmm::clear(U1);
				gmm::copy(UM_transp, U1);
				cout << endl<<"... time to solve U1: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A2: 
			{	gmm::resize(U2, dof_tot);	gmm::clear(U2);
				gmm::copy(UM_transp, U2);	
				cout << endl<<"... time to solve U2: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A3: 
			{	gmm::resize(U3, dof_tot);	gmm::clear(U3);
				gmm::copy(UM_transp, U3);	
				cout << endl<<"... time to solve U3: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
		 }
				break;} 

		case DUAL:{
		switch(ASSUMPTION){

			case A0: 
			{	gmm::resize(Z0, dof_tot);	gmm::clear(Z0);
				gmm::copy(UM_transp, Z0);	
				cout << endl<<"... time to solve Z0: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A1: 
			{	gmm::resize(Z1, dof_tot);	gmm::clear(Z1);
				gmm::copy(UM_transp, Z1);	
				cout << endl<<"... time to solve Z1: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A2: 
			{	gmm::resize(Z2, dof_tot);	gmm::clear(Z2);
				gmm::copy(UM_transp, Z2);	
				cout << endl<<"... time to solve Z2: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A3: 
			{	gmm::resize(Z3, dof_tot);	gmm::clear(Z3);
				gmm::copy(UM_transp, Z3);	
				cout << endl<<"... time to solve Z3: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
		 }
				break;} 

		}
	return true;

	};
	//! Export problem U1
	void transport3d1d::export_model(const size_type ASSUMPTION,const  size_type VERSION, const string & suff ){

	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	// define dofs
	size_type dof_t;
	size_type dof_v;

	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		dof_t  = dof_transp.Ct_Omega();	}
	else{
		dof_t  = dof_transp.Ct();	}

 	if(ASSUMPTION==A0){
		dof_v  = dof_transp.Ct_Sigma();	}
	else{
		dof_v  = dof_transp.Cv();	}

	size_type dof_tot= dof_t+dof_v;
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_t); 

	// Array of unknown dof of the network velocity
	vector_type Cv(dof_v); 

	//Copy solution


	switch(VERSION){

		case PRIMAL:{
		switch(ASSUMPTION){

			case A0: 
			{	gmm::copy(gmm::sub_vector(U0, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U0, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A1: 
			{	gmm::copy(gmm::sub_vector(U1, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U1, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A2: 
			{	gmm::copy(gmm::sub_vector(U2, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U2, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A3: 
			{	gmm::copy(gmm::sub_vector(U3, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U3, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
		}
				break;} 

		case DUAL:{
		switch(ASSUMPTION){

			case A0: 
			{	gmm::copy(gmm::sub_vector(Z0, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z0, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A1: 
			{	gmm::copy(gmm::sub_vector(Z1, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z1, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A2: 
			{	gmm::copy(gmm::sub_vector(Z2, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z2, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A3: 
			{	gmm::copy(gmm::sub_vector(Z3, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z3, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
		}
				break;} 

	}

	string assumption, version;	
	switch(VERSION){
		case PRIMAL:
			version="U";
			break;
		case DUAL:
			version="Z";
			break;
	}
	
	switch(ASSUMPTION){
		case A0:
			assumption="0";
			break;
		case A1:
			assumption="1";
			break;
		case A2:
			assumption="2";
			break;
		case A3:
			assumption="3";
			break;
	}
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting reduced solution in tissue: "<<version<<assumption << endl;
	#endif
	if(ASSUMPTION==A1||ASSUMPTION==A0){
		vtk_export exp_Ct(descr_transp.OUTPUT+version+"t"+assumption+suff+".vtk");
		exp_Ct.exporting(mf_Ct_Omega);
		exp_Ct.write_mesh();
		exp_Ct.write_point_data(mf_Ct_Omega, Ct, version+"t"+assumption);
	}else{
		vtk_export exp_Ct(descr_transp.OUTPUT+version+"t"+assumption+suff+".vtk");
		exp_Ct.exporting(mf_Ct);
		exp_Ct.write_mesh();
		exp_Ct.write_point_data(mf_Ct, Ct, version+"t"+assumption);
	}


	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting reduced solution in vessel: "<<version<<assumption << endl;
	#endif

	if(ASSUMPTION==A0){
		vtk_export exp_Cv(descr_transp.OUTPUT+version+"v"+assumption+suff+".vtk");
		exp_Cv.exporting(mf_Ct_Sigma);
		exp_Cv.write_mesh();
		exp_Cv.write_point_data(mf_Ct_Sigma, Cv, version+"v"+assumption);
		
	}else{
		vtk_export exp_Cv(descr_transp.OUTPUT+version+"v"+assumption+suff+".vtk");
		exp_Cv.exporting(mf_Cv);
		exp_Cv.write_mesh();
		exp_Cv.write_point_data(mf_Cv, Cv, version+"v"+assumption);
	}

	if(ASSUMPTION==A0){
		vector_type u0(dof_transp.Ct_Sigma()); 	
		vector_type u0_1D(dof_transp.Cv());
		vector_type u0_3D(dof_transp.Ct());

		gmm::copy(gmm::sub_vector(U0, 
				  gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Ct_Sigma())), u0);

		getfem::interpolation(mf_Ct_Sigma, mf_Ct, u0, u0_3D);

		gmm::mult(MBARBAR, u0_3D,u0_1D);		
		vtk_export exp_Cv2(descr_transp.OUTPUT+"Uv0_1D.vtk");
		exp_Cv2.exporting(mf_Cv);
		exp_Cv2.write_mesh();
		exp_Cv2.write_point_data(mf_Cv, u0_1D, "Uv0_1D.vtk");
	}


	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl<<endl; 
	#endif
};


	/////////////////////////////////////////////
	// Compute model errors

	//! Compute model error of A1	
	void transport3d1d::compute_model_error_1(void){
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A1 ..." << endl;
	#endif
	// Storing error model from bilinear form 
	vector_type d1(dof_transp.Cv()); gmm::clear(d1);
	cout<<"WARNING: d1(u,z)=0!"<<endl;

	// Storing error model from linear form 
	vector_type l1_t(dof_transp.Cv()); gmm::clear(l1_t);
	vector_type l1_g(dof_transp.Cv()); gmm::clear(l1_g);
	//cout<<"WARNING: l1(u,z)=0 if g=const!!"<<endl;


	// Solution of dual 1D problem
	vector_type z1(dof_transp.Cv()); gmm::clear(z1);
	gmm::copy(gmm::sub_vector(Z1, 
		gmm::sub_interval(0, dof_transp.Cv())), z1);

	// Solution of primal 3D problem	
	vector_type u1(dof_transp.Ct_Omega()); gmm::clear(u1);

	gmm::copy(gmm::sub_vector(U0, 
				  gmm::sub_interval(0, dof_transp.Ct_Omega())), u1);

	gmm::add(gmm::scaled(gmm::sub_vector(U1, 
				  gmm::sub_interval(0, dof_transp.Ct_Omega())),-1), u1);


	// Build l1_t, due to the difference of ut0 and ut1 on Gamma

	//Perimeter = 2*pi*R*k 
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);

/*	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Mvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mvv);
	
	getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Perimeter);
	gmm::mult(MBAR, Mvv, Bvt); */
	
	// l1_t = z1 *(Bvt*u1)		
	gmm::mult(BVT_Omega, u1, l1_t);   
	gmm::vscale(l1_t, z1);

	// Build l1_g, due to the difference of g and his average on a cross section
	
//  do nothing, l1_g is zero

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting d1 and l1 ..." << endl;
	#endif
	vtk_export exp_d(descr_transp.OUTPUT+"d1.vtk");
	exp_d.exporting(mf_Cv);
	exp_d.write_mesh();
	exp_d.write_point_data(mf_Cv, d1, "d1");

	vtk_export exp_lt(descr_transp.OUTPUT+"l1.vtk");
	exp_lt.exporting(mf_Cv);
	exp_lt.write_mesh();
	exp_lt.write_point_data(mf_Cv, l1_t, "l1");

	vtk_export exp_lg(descr_transp.OUTPUT+"l1_g.vtk");
	exp_lg.exporting(mf_Cv);
	exp_lg.write_mesh();
	exp_lg.write_point_data(mf_Cv, l1_g, "l1_g");
	
	};

	//! Compute model error of A2
	void transport3d1d::compute_model_error_2(void){

	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A2 ..." << endl;
	#endif
	// Storing error model from bilinear form 
	vector_type d2(dof_transp.Ct()); gmm::clear(d2);
	// Storing error model from linear form 
	vector_type l2(dof_transp.Ct()); gmm::clear(l2);
	// Storing error model from linear form 
	vector_type l2_v(dof_transp.Ct()); gmm::clear(l2_v);

	// Solution of primal 3D problem
	vector_type u2(dof_transp.Ct()); gmm::clear(u2);
	// Solution of dual 3D problem
	vector_type z2(dof_transp.Ct()); gmm::clear(z2);
	// Source for primal 3D problem
	scalar_type f = PARAM.real_value("F", "Value of source for tissue reduced problem");	
	vector_type F(dof.coeft(),f);

	gmm::copy(gmm::sub_vector(U2, 
		gmm::sub_interval(0, dof_transp.Ct())), u2);
	gmm::copy(gmm::sub_vector(Z2, 
		gmm::sub_interval(0, dof_transp.Ct())), z2);

	//assemble estimator d2
	sparse_matrix_type At(dof_transp.Ct(), dof_transp.Ct());gmm::clear(At);
	getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At, mimt, mf_Ct, descr_transp.SIGMA);
	gmm::mult(At, u2, d2);   // M1 * V2 --> V1
	gmm::vscale(d2, z2);
	gmm::scale(d2, -1);

	//assemble estimator l2
	getfem::asm_source_term(l2, mimt, mf_Ct, mf_coeft, F, descr_transp.SIGMA);
	gmm::vscale(l2, z2);
	gmm::scale(l2, -1);

	//assemble estimator l2_1
	vector_type uv(dof_transp.Cv()); gmm::clear(uv);
	gmm::copy(gmm::sub_vector(U1, 
				  gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), uv);

	gmm::add(gmm::scaled(gmm::sub_vector(U2, 
				  gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())),-1), uv);

	//Perimeter = 2*pi*R*k 
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);

	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Mtv(dof_transp.Ct(), dof_transp.Cv());gmm::clear(Mtv);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Mvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mvv);
	
	getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Perimeter);
	gmm::mult(gmm::transposed(MBAR), Mvv, Mtv); // M1 * M2 ---> M3	
	gmm::mult(Mtv, uv, l2_v);
	gmm::vscale(l2_v, z2);
	
	//export d2 and l2

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting d2 and l2 ..." << endl;
	#endif
	vtk_export exp_d(descr_transp.OUTPUT+"d2.vtk");
	exp_d.exporting(mf_Ct);
	exp_d.write_mesh();
	exp_d.write_point_data(mf_Ct, d2, "d2");

	vtk_export exp_l(descr_transp.OUTPUT+"l2.vtk");
	exp_l.exporting(mf_Ct);
	exp_l.write_mesh();
	exp_l.write_point_data(mf_Ct, l2, "l2");

	vtk_export exp_lv(descr_transp.OUTPUT+"l2_v.vtk");
	exp_lv.exporting(mf_Ct);
	exp_lv.write_mesh();
	exp_lv.write_point_data(mf_Ct, l2_v, "l2_v");

	};

	//! Compute model error of A3	
	void transport3d1d::compute_model_error_3(void){
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A3 ..." << endl;
	#endif

	// Storing error model from bilinear form 
	vector_type d3(dof_transp.Ct()); gmm::clear(d3);

	//Storing error model from linear form 
	vector_type l3(dof_transp.Ct()); gmm::clear(l3);
	//cout<<"WARNING: l3(u,z)=0!!"<<endl;

	// Solution of primal 3D problem
	vector_type u3(dof_transp.Ct()); gmm::clear(u3);
	// Solution of dual 3D problem
	vector_type z3(dof_transp.Ct()); gmm::clear(z3);
	gmm::copy(gmm::sub_vector(U3, 
		gmm::sub_interval(0, dof_transp.Ct())), u3);
	gmm::copy(gmm::sub_vector(Z3, 
		gmm::sub_interval(0, dof_transp.Ct())), z3);

	//Assemble estimator d3
	sparse_matrix_type Mtt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Mtt);
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	vector_type Kt(dof.coeft(),k);
	getfem::asm_mass_matrix_param (Mtt, mimt, mf_Ct,mf_coeft, Kt, descr_transp.GAMMA);
	gmm::scale(Mtt,-1);
	gmm::add(BTT,Mtt);
	gmm::mult(Mtt, u3, d3);   // M1 * V2 --> V1
	gmm::vscale(d3, z3);

	//assemble estimator l3

	//assemble estimator l2_1
	vector_type uv(dof_transp.Cv()); gmm::clear(uv);
	gmm::copy(gmm::sub_vector(U2, 
				  gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), uv);

	gmm::add(gmm::scaled(gmm::sub_vector(U3, 
				  gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())),-1), uv);

	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);

	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Mtv(dof_transp.Ct(), dof_transp.Cv());gmm::clear(Mtv);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Mvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mvv);
	
	getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Perimeter);
	gmm::mult(gmm::transposed(MBAR), Mvv, Mtv); // M1 * M2 ---> M3	
	gmm::mult(Mtv, uv, l3);
	gmm::vscale(l3, z3);

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting d3 and l3 ..." << endl;
	#endif
	vtk_export exp_d(descr_transp.OUTPUT+"d3.vtk");
	exp_d.exporting(mf_Ct);
	exp_d.write_mesh();
	exp_d.write_point_data(mf_Ct, d3, "d3");

	vtk_export exp_l(descr_transp.OUTPUT+"l3.vtk");
	exp_l.exporting(mf_Ct);
	exp_l.write_mesh();
	exp_l.write_point_data(mf_Ct, l3, "l3");


	};


	/////////////////////////////////////////////
	//User interface

	//User function to compute all the the steps (reduced primal and Z0,Z1,Z2)
	void transport3d1d::model_error(int argc, char *argv[]){

	double timeA0=0; 
	double timeA1=0;
	double timeA2=0;
	double timeA3p=0;
	double timeA3d=0;
	double time_model=0;


 	double time = gmm::uclock_sec();

		//initialize 
		init_model(argc, argv);

/*	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"Z0: solve reference problem:"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	
		//assemble        
		assembly_model(A0);   
		assembly_dual_model(A0);    
		//solve     
		if (!solve_model(A0,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A0,DUAL);

	timeA0= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();
*/	
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"Z1: solve problem for first assumption:"<<endl;
	cout <<"    U(s,r,theta)=U(s)"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	 
		//assemble        
		assembly_model(A1); 
		assembly_dual_model(A1);    
		//solve     
		if (!solve_model(A1,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,DUAL);

	timeA1= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"Z2: solve problem for second assumption:"<<endl;
	cout <<"    Omega+ == Omega"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

		//assemble        
		assembly_model(A2); 
		assembly_dual_model(A2);    
		//solve     
		if (!solve_model(A2,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,DUAL);


	timeA2= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble      
		assembly_model(A3);    
		//solve     
		if (!solve_model(A3,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,PRIMAL);

	timeA3p= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"Z3: solve problem for third assumption:"<<endl;
	cout <<"    Neglect fluctuations"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

		//assemble        
		assembly_dual_model(A3);    
		//solve     
		if (!solve_model(A3,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,DUAL);


	timeA3d= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D1 AND L1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		//compute_model_error_1();

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D2 AND L2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		//compute_model_error_2();

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D3 AND L3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		//compute_model_error_3();

	time_model= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	cout << "========================================================="<<endl;
	cout << endl<<"... time to solve Z0: " << timeA0 << " seconds\n";
	cout << endl<<"... time to solve Z1: " << timeA1 << " seconds\n";
	cout << endl<<"... time to solve Z2: " << timeA2 << " seconds\n";
	cout << endl<<"... time to solve Z3: " << timeA3d << " seconds\n";
	cout << endl<<"... time to solve U3: " << timeA3p << " seconds\n";
	cout << endl<<"... time to solve error estimators: " << time_model << " seconds\n";

};



	//User function to compute all the the steps
	void transport3d1d::model_error_old(int argc, char *argv[]){

	double time = gmm::uclock_sec();
	double timeA0=0; 
	double timeA1=0;
	double timeA2=0;
	double timeA3=0;

 	
		//initialize 
		init_model(argc, argv);

	
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A0: solve reference problem:"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  


	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U0 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A0);    
		//solve     
		if (!solve_model(A0,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A0,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z0 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A0);    
		//solve     
		if (!solve_model(A0,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A0,DUAL);

	timeA0= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();
	
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A1: solve problem for first assumption:"<<endl;
	cout <<"    U(s,r,theta)=U(s)"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A1);    
		//solve     
		if (!solve_model(A1,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A1);    
		//solve     
		if (!solve_model(A1,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D1 AND L1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_1();

	timeA1= time-gmm::uclock_sec();
	time= gmm::uclock_sec();

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A2: solve problem for second assumption:"<<endl;
	cout <<"    Omega+ == Omega"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A2);    
		//solve     
		if (!solve_model(A2,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A2);    
		//solve     
		if (!solve_model(A2,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D2 AND L2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_2();

	timeA2= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A3: solve problem for third assumption:"<<endl;
	cout <<"    Neglect fluctuations"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A3);    
		//solve     
		if (!solve_model(A3,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A3);    
		//solve     
		if (!solve_model(A3,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D3 AND L3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_3();

	timeA3= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	cout << "========================================================="<<endl;
	cout << endl<<"... time to solve A0: " << timeA0 << " seconds\n";
	cout << endl<<"... time to solve A1: " << timeA1 << " seconds\n";
	cout << endl<<"... time to solve A2: " << timeA2 << " seconds\n";
	cout << endl<<"... time to solve A3: " << timeA3 << " seconds\n";

};

	//User method to compute model error for A1
	void transport3d1d::model_error_A0(int argc, char *argv[]){// A1: U=U(s)
	//Primal	
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A0: solve reference problem:"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

		//initialize 
		init_model(argc, argv);


	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U0 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A0);    
		//solve     
		if (!solve_model(A0,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A0,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z0 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A0);    
		//solve     
		if (!solve_model(A0,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A0,DUAL);




};

	//User method to compute model error for A1
	void transport3d1d::model_error_A1(int argc, char *argv[]){// A1: U=U(s)
	//Primal	
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A1: solve problem for first assumption:"<<endl;
	cout <<"    U(s,r,theta)=U(s)"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

		//initialize 
		init_model(argc, argv);


	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A1);    
		//solve     
		if (!solve_model(A1,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A1);    
		//solve     
		if (!solve_model(A1,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D1 AND L1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_1();


};


	//User method to compute model error for A2
	void transport3d1d::model_error_A2(int argc, char *argv[]){//A2: Omega+ == Omega
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A2: solve problem for second assumption:"<<endl;
	cout <<"    Omega+ == Omega"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;

		//initialize 
		init_model(argc, argv);

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A2);    
		//solve     
		if (!solve_model(A2,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A2);    
		//solve     
		if (!solve_model(A2,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D2 AND L2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_2();
		

	};


	//User method to compute model error for A3
	void transport3d1d::model_error_A3(int argc, char *argv[]){//A3: neglect fluctuations

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A3: solve problem for third assumption:"<<endl;
	cout <<"    Neglect fluctuations"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;

		//initialize 
		init_model(argc, argv);
	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A3);    
		//solve     
		if (!solve_model(A3,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A3);    
		//solve     
		if (!solve_model(A3,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D3 AND L3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_3();

	};



} // end of namespace
