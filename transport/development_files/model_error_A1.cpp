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

const int A1=1;
const int A2=2;
const int A3=3;

const int PRIMAL=1;
const int DUAL=2;


	//Aux methods for model error
	//! Initialize problem U1
	void init_U1(int argc, char *argv[]);
	//! Assemble problem U1
	bool assembly_U1(void);
	//! Solve problem U1
	bool solve_U1(void);
	//! Export problem U1
	void export_U1(const string & suff = "");
	//! Initialize problem Z1
	void init_Z1(int argc, char *argv[]);
	//! Assemble problem Z1
	bool assembly_Z1(void);
	//! Solve problem Z1
	bool solve_Z1(void);
	//! Export problem Z1
	void export_Z1(const string & suff = "");
	//! Compute model error of A1	
	void compute_error_1(void);

	//! Initialize problem U2
	void init_U2(int argc, char *argv[]);
	//! Assemble problem U2
	bool assembly_U2(void);
	//! Solve problem U2
	bool solve_U2(void);
	//! Export problem U2
	void export_U2(const string & suff = "");
	//! Initialize problem Z2
	void init_Z2(int argc, char *argv[]);
	//! Assemble problem Z2
	bool assembly_Z2(void);
	//! Solve problem Z2
	bool solve_Z2(void);
	//! Export problem Z2
	void export_Z2(const string & suff = "");
	//! Compute model error of A2
	void compute_error_2(void);

	//! Initialize problem U3
	void init_U3(int argc, char *argv[]);
	//! Assemble problem U3
	bool assembly_U3(void);
	//! Solve problem U3
	bool solve_U3(void);
	//! SoExportlve problem U3
	void export_U3(const string & suff = "");
	//! Initialize problem Z3
	void init_Z3(int argc, char *argv[]);
	//! Assemble problem Z3
	bool assembly_Z3(void);
	//! Solve problem Z3
	bool solve_Z3(void);
	//! Export problem Z3
	void export_Z3(const string & suff = "");
	//! Compute model error of A3
	void compute_error_3(void);


	//! Initialize problem U1
	void transport3d1d::init_U1(int argc, char *argv[]){

	PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	import_data_transp();
	//2. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh_transp();
	//3. Set finite elements and integration methods
	set_im_and_fem_U1();
 	//4. Build problem parameters
	build_param_transp();
	//5. Build the list of tissue boundary data
	build_tissue_boundary_transp();
	//6. Build the list of tissue boundary (and junction) data
 	build_vessel_boundary_transp();
};
	//3. Set finite elements and integration methods: mf_t should be defined only on Omega_plus
	void transport3d1d::set_im_and_fem_U1(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	
	
	pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
	pfem pf_Cv = fem_descriptor(descr_transp.FEM_TYPEV_C); 

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for tissue ..." << endl;
	#endif
		


	mf_Ct.set_finite_element(mesht.convex_index(), pf_Ct);
	mf_Ct_Omega.set_finite_element(mesht.region(descr_transp.OMEGA).index(), pf_Ct); //.region(descr_transp.OMEGA)
	mf_Ct_Sigma.set_finite_element(mesht.region(descr_transp.SIGMA).index(), pf_Ct); //.region(descr_transp.SIGMA)

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


}



;
	//! Assemble problem U1
	bool transport3d1d::assembly_U1(void){

	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM_transp, UM_transp, FM_transp ..." << endl;
	#endif
	gmm::resize(AM_transp, dof_transp.Ct_Omega()+dof_transp.Cv(), 
			       dof_transp.Ct_Omega()+dof_transp.Cv());	gmm::clear(AM_transp);
	gmm::resize(UM_transp, dof_transp.Ct_Omega()+dof_transp.Cv());	gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_transp.Ct_Omega()+dof_transp.Cv()); 	gmm::clear(FM_transp);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM_transp ..." << endl;
	#endif
	// Mass(time derivative)  matrix for the interstitial problem
	sparse_matrix_type At(dof_transp.Ct_Omega(), dof_transp.Ct_Omega());gmm::clear(At);
	// Mass (time derivative)  matrix for the network problem
	sparse_matrix_type Av(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Av);	

	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Mtt(dof_transp.Ct_Omega(), dof_transp.Ct_Omega());gmm::clear(Mtt);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Mtv(dof_transp.Ct_Omega(), dof_transp.Cv());gmm::clear(Mtv);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Mvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mvv);

	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof_transp.Ct_Omega(), dof_transp.Ct());gmm::clear(Btt);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof_transp.Ct_Omega(), dof_transp.Cv());gmm::clear(Btv);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_transp.Cv(), dof_transp.Ct_Omega());gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bvv);
	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar(dof_transp.Cv(), dof_transp.Ct_Omega());gmm::clear(Mbar);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof_transp.Cv(), dof_transp.Ct_Omega());gmm::clear(Mlin);
	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbarbar(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mbarbar);
	
	// Aux tissue source vector
	vector_type Ft(dof_transp.Ct_Omega()); gmm::clear(Ft);
	// Aux vessels source vector
	vector_type Fv(dof_transp.Cv()); gmm::clear(Fv);

	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix At ..." << endl;
	#endif	
	// Assemble At, stiffness matrix for laplacian in tissue
	getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct_Omega); 
	
	gmm::add(At,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct_Omega()), 
					gmm::sub_interval(0, dof_transp.Ct_Omega()))); 


	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix Av ..." << endl;
	#endif	
	// Assemble Av, stiffness matrix for laplacian in vessel
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); //Area = pi*R^2

	getfem::asm_stiffness_matrix_for_laplacian(Av,mimv,mf_Cv,mf_coefv,Area); 
	
	gmm::add(Av,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv()))); 



	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	vector_type Kv(dof.coefv(),k);
	vector_type Kt(dof.coeft(),k);

	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);
	bool NEWFORM = 1;//PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");

	asm_exchange_aux_mat_transp(Mbar, Mlin, 
			mimv, mf_Ct_Omega, mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);	
	asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, Mbar, Mlin, Perimeter, NEWFORM);

	// Copying Bvt

	gmm::add(scaled(Bvt,-1),								
			  gmm::sub_matrix(AM_transp, 
			  		gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv()),
					gmm::sub_interval(0, dof_transp.Ct_Omega())));
	// Copying Bvv

	gmm::add(Bvv,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv()))); 
	

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling mass matrices on Gamma..." << endl;
	#endif
	// Assemble Mtt, mass matrix on tissue
	
	// Assemble Mtv, exchange matrix on gamma
	getfem::asm_mass_matrix_param (Mtt, mimt, mf_Ct_Omega,mf_coeft, Kt, descr_transp.GAMMA);
	gmm::vscale(Area, Kv);

	getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Kv);

	gmm::mult(gmm::transposed(Mbar), Mvv, Mtv); // M1 * M2 ---> M3	
	
	// Copying Mtt
	gmm::add(Mtt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct_Omega()), 
					gmm::sub_interval(0, dof_transp.Ct_Omega()))); 
	// Copying Mtv

	gmm::add(scaled(Mtv,-1),
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct_Omega()),
					gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv()))); 


	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling source terms ..." << endl;
	#endif

	// Assemble F: source term in tissue
	scalar_type f = PARAM.real_value("F", "Value of source for tissue reduced problem");	
	vector_type F(dof.coeft(),f);
	getfem::asm_source_term(Ft, mimt, mf_Ct_Omega, mf_coeft, F);
	
	// Assemble G: source term in vessel
	//! /todo Extend to variable source term: use muParser and build cross section average (operator Mbarbar)
	bool G_CONSTANT = PARAM.int_value("G_CONSTANT", "flag for using constant source term g in vessel");
	if(G_CONSTANT){
	scalar_type g = PARAM.real_value("G", "Value of source for vessel reduced problem");
	vector_type G(dof.coefv(),g);	
	gmm::vscale(Area, G); //G=G*pi*R^2
	getfem::asm_source_term(Fv, mimv, mf_Cv, mf_coefv, G);
	}
	else{
	scalar_type g = PARAM.real_value("G", "Value of source for vessel reduced problem");
	vector_type G(dof.coeft(),g);	
	vector_type Areat(dof.coeft());
	getfem::interpolation(mf_coefv, mf_coeft, Area, Areat);
	gmm::vscale(Areat, G); //G=G*pi*R^2

	// Aux tissue source vector
	vector_type Ftt(dof_transp.Ct()); gmm::clear(Ftt);
	getfem::asm_source_term(Ftt, mimt, mf_Ct, mf_coeft, G);

	//Build average matrix Mbarbar
	asm_exchange_aux_mat_bar_bar(Mbarbar,mimt, mf_Ct, mf_Cv,mf_coefv, param.R(), descr.NInt, descr_transp.NIntA, nb_branches);

	gmm::mult(Mbarbar, Ftt, Fv);
	}

	gmm::add(Ft, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(0,dof_transp.Ct_Omega())));

	gmm::add(Fv, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv())));


	// De-allocate memory
	gmm::clear(At);	   gmm::clear(Av); 
	gmm::clear(Mtt);   gmm::clear(Mtv); 	    
	gmm::clear(Mbar);  gmm::clear(Mlin);
	gmm::clear(Btt);   gmm::clear(Btv);
	gmm::clear(Bvt);   gmm::clear(Bvv);

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling boundary conditions ..." << endl;
	#endif	
	gmm::resize(AM_temp, dof_transp.Ct_Omega()+dof_transp.Cv(),
			     dof_transp.Ct_Omega()+dof_transp.Cv());	gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_transp.Ct_Omega()+dof_transp.Cv()); 	gmm::clear(FM_temp);
	gmm::copy(AM_transp,AM_temp);
	gmm::copy(FM_transp,FM_temp);

	vector_type Dirichlet_null(dof_transp.Ct_Omega(),0);
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
		getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct_Omega, BCt_transp[bc].rg, Dirichlet_null);	//Homogeneous Dirichlet on all the faces of tissue			
	} 
	
	//On vessel, we have homogeneous neumann conditions: we simply don't implement the boundary term due to diffusion.

};
	//! Solve problem U1
	bool transport3d1d::solve_U1(void){

  	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	double time = gmm::uclock_sec();
	gmm::clear(UM_transp);
	gmm::clean(AM_temp, 1E-12);
	gmm::clean(FM_temp, 1E-12);
	// Solve the system on AM_temp, UM_transp, FM_temp
	bool solved = solver_transp();
	if (!solved) return false;
	
	gmm::resize(U1, dof_transp.Ct_Omega()+dof_transp.Cv());	gmm::clear(U1);
	gmm::copy(UM_transp, U1);
	//export solution
	#ifdef M3D1D_VERBOSE_
	std::cout<<"solved! going to export..."<<std::endl;
	#endif	
	
	cout << endl<<"... time to solve U1: " << gmm::uclock_sec() - time << " seconds\n";				
	
	return true;

};
	//! Export problem U1
	void transport3d1d::export_U1(const string & suff ){

	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_transp.Ct_Omega()); 

	// Array of unknown dof of the network velocity
	vector_type Cv(dof_transp.Cv()); 

	//Copy solution
	gmm::copy(gmm::sub_vector(U1, 
		gmm::sub_interval(0, dof_transp.Ct_Omega())), Ct);
	gmm::copy(gmm::sub_vector(U1, 
		gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);

	

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ct ..." << endl;
	#endif
	vtk_export exp_Ct(descr_transp.OUTPUT+"Ut1.vtk");
	exp_Ct.exporting(mf_Ct_Omega);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct_Omega, Ct, "Ut1");



	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+"Uv1.vtk");
	exp_Cv.exporting(mf_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_Cv, Cv, "Uv1");

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif
};



	//! Initialize problem Z1
	void transport3d1d::init_Z1(int argc, char *argv[]){

	// Do nothing! You should already have implemented init_U1 although
};




	//! Assemble problem Z1
	bool transport3d1d::assembly_Z1(void){

	gmm::clear(UM_transp);

	// Aux tissue source vector
	vector_type Jt(dof_transp.Ct_Omega()); gmm::clear(Jt);
	// Aux vessels source vector
	vector_type Jv(dof_transp.Cv()); gmm::clear(Jv);

	// Aux tissue source vector
	vector_type Ft(dof.coeft(),1);
	// Aux vessels source vector
	vector_type Fv(dof.coefv(),1);
 


	// prepare AM and Fm temp
	gmm::resize(AM_temp, dof_transp.Ct_Omega()+dof_transp.Cv(),
			     dof_transp.Ct_Omega()+dof_transp.Cv());	gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_transp.Ct_Omega()+dof_transp.Cv()); 	gmm::clear(FM_temp);
	gmm::copy(AM_transp,AM_temp);
	gmm::copy(FM_transp,FM_temp);

	gmm::scale(gmm::sub_matrix(AM_temp, 
					gmm::sub_interval(0, dof_transp.Ct_Omega()),
					gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv()))
		 ,0);
	gmm::scale(gmm::sub_matrix(AM_temp, 
			  		gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv()),
					gmm::sub_interval(0, dof_transp.Ct_Omega()))
		 ,0);

		#ifdef M3D1D_VERBOSE_
	cout << "  Assembling source terms ..." << endl;
	#endif
	// Assemble F: source term in tissue
	getfem::asm_source_term(Jt,mimt, mf_Ct_Omega, mf_coeft, Ft);

	// Assemble G: source term in vessel
	getfem::asm_source_term(Jv,mimv, mf_Cv, mf_coefv, Fv ); 

	//Add to FM_temp

	gmm::add(Jt, 
			gmm::sub_vector(FM_temp,
					gmm::sub_interval(0,dof_transp.Ct_Omega())));

	gmm::add(Jv, 
			gmm::sub_vector(FM_temp,
					gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv())));

/*
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

	sparse_matrix_type JJv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(JJv);
		
	generic_assembly L2normv;
	L2normv.push_mi(mimv);
	L2normv.push_mf(mf_Cv);
	L2normv.push_mat(JJv);
	L2normv.set("V(#1,#1)+=comp(Base(#1).Base(#1))(:,:)");
	L2normv.assembly();
	
	//gmm::add(gmm::mat_trace(JJv),Jv);
	for (int i=0; i>dof_transp.Cv(); i++){
		Jv[i]=JJv(i,i);
		Jv[i]=sqrt(Jv[i]);
		
	};
*/	





	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling boundary conditions ..." << endl;
	#endif	

	vector_type Dirichlet_null(dof_transp.Ct_Omega(),0);
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
		getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct_Omega, BCt_transp[bc].rg, Dirichlet_null);	//Homogeneous Dirichlet on all the faces of tissue			
	} 
};
	//! Solve problem Z1
	bool transport3d1d::solve_Z1(void){

 	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	double time = gmm::uclock_sec();
	gmm::clear(UM_transp);
	gmm::clean(AM_temp, 1E-12);
	gmm::clean(FM_temp, 1E-12);
	// Solve the system on AM_temp, UM_transp, FM_temp
	bool solved = solver_transp();
	if (!solved) return false;
	
	gmm::resize(Z1, dof_transp.Ct_Omega()+dof_transp.Cv());	gmm::clear(Z1);
	gmm::copy(UM_transp, Z1);
	//export solution
	#ifdef M3D1D_VERBOSE_
	std::cout<<"solved! going to export..."<<std::endl;
	#endif	
	
	cout << endl<<"... time to solve Z1: " << gmm::uclock_sec() - time << " seconds\n";				
	
	return true;

};
	//! Export problem Z1
	void transport3d1d::export_Z1(const string & suff ){

	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_transp.Ct_Omega()); 

	// Array of unknown dof of the network velocity
	vector_type Cv(dof_transp.Cv()); 

	//Copy solution
	gmm::copy(gmm::sub_vector(Z1, 
		gmm::sub_interval(0, dof_transp.Ct_Omega())), Ct);
	gmm::copy(gmm::sub_vector(Z1, 
		gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Cv())), Cv);

	

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ct ..." << endl;
	#endif
	vtk_export exp_Ct(descr_transp.OUTPUT+"Zt1.vtk");
	exp_Ct.exporting(mf_Ct_Omega);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct_Omega, Ct, "Zt1");



	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+"Zv1.vtk");
	exp_Cv.exporting(mf_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_Cv, Cv, "Zv1");

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif
};

	//! Compute model error of A1	
	void transport3d1d::compute_error_1(void){
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A1 ..." << endl;
	#endif
	// Storing error model from bilinear form 
	vector_type d1(dof_transp.Cv()); gmm::clear(d1);
	cout<<"WARNING: d1(u,z)=0!"<<endl;

	// Storing error model from linear form 
	vector_type l1(dof_transp.Cv()); gmm::clear(l1);
	cout<<"WARNING: l1(u,z)=0 if g=const!!"<<endl;

// do nothing: estimators are zero!


	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting d1 and l1 ..." << endl;
	#endif
	vtk_export exp_d(descr_transp.OUTPUT+"d1.vtk");
	exp_d.exporting(mf_Cv);
	exp_d.write_mesh();
	exp_d.write_point_data(mf_Cv, d1, "d1");

	vtk_export exp_l(descr_transp.OUTPUT+"l1.vtk");
	exp_l.exporting(mf_Cv);
	exp_l.write_mesh();
	exp_l.write_point_data(mf_Cv, l1, "l1");
	
};


 
 } // end of namespace
