/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017

                  
                Copyright (C) 2016 Stefano Brambilla
  ======================================================================
 "Computational models for nanoparticle transport in the vascular system"
    	       Master thesis in Mathematical Engineering
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
         Copyright (C) 2017 Annagiulia Tiozzo, Federica Laurino
======================================================================*/
/*! 
  @file   transport3d1d.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
	  Annagiulia Tiozzo <annagiulia92t@gmail.com>
	  Federica Laurino <federica.laurino@polimi.it>
  @date   April 2017.
  @brief  Definition of the main class for the 3D/1D coupled transport problem.
 */
 
 #include <transport3d1d.hpp>
 #include <assembling1d_transp_nano.hpp>

 namespace getfem {
 
 void transport3d1d::init (int argc, char *argv[]) 
 {
 	std::cout << "initialize transport problem..."<<std::endl<<std::endl;
 
 	import_data();
 	build_mesh();
 	set_im_and_fem();
 	build_param();
 	build_vessel_boundary();
	build_tissue_boundary();
 
 }; // end of init


 // Aux methods for init
	
//! Import algorithm specifications
void transport3d1d::import_data(void)
{
	std::cout<<"init part 1: import data!......" <<std::endl;
	#ifdef M3D1D_VERBOSE_
	cout << "Importing descriptors for tissue and vessel problems ..." << endl;
	#endif
	descr_transp.import(PARAM);
	#ifdef M3D1D_VERBOSE_
	cout << descr_transp;
	#endif
	 
 
};
	
	 
	
//! Import mesh for tissue (3D) and vessel (1D)  
void transport3d1d::build_mesh(void)
{
	//but, in order to have the boundary conditions for the nodes
	//we need to build again the 1D mesh from another pts file
	mesht.clear();
		bool test = 0;
	test = PARAM.int_value("TEST_GEOMETRY");
	if(test==0){
		#ifdef M3D1D_VERBOSE_
		cout << "Importing the 3D mesh for the tissue ...  "   << endl;
		#endif
		 import_msh_file(descr.MESH_FILET, mesht);
	}else{
		#ifdef M3D1D_VERBOSE_
		cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
		#endif
		string st("GT='" + PARAM.string_value("GT_T") + "'; " +
					   "NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
					   "ORG=" + PARAM.string_value("ORG_T") + "; " +  
					   "SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
					   "NOISED=" + PARAM.string_value("NOISED_T")); 
		cout << "mesht description: " << st << endl;
		regular_mesh(mesht, st);
	}
	
	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel (transport problem)... "   << endl;
	#endif
	std::ifstream ifs(descr_transp.MESH_FILEV);
	meshv.clear();

	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr_transp.MESH_FILEV);
	import_pts_file(ifs, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV);
	nb_branches = nb_vertices.size();
	ifs.close();
	
	
};

//! Set finite elements methods and integration methods 
void transport3d1d::set_im_and_fem(void)
{
	std::cout<<"init part 2: set fem methods!......" <<std::endl;
	

	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	
	
	pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
	pfem pf_Cv = fem_descriptor(descr_transp.FEM_TYPEV_C);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for tissue ..." << endl;
	#endif
		

	mf_Ct.set_finite_element(mesht.convex_index(), pf_Ct);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif

	mf_Cv.set_finite_element(meshv.convex_index(), pf_Cv);

	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif

	dof_transp.set(mf_Ct, mf_Cv);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof_transp;
	#endif
	

	//Define once again the finite element methods and the integration method on the network
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
	
	
//! Build problem parameters
void transport3d1d::build_param(void)
{
	std::cout<<"init part 3: build dimensionless parameters!" <<std::endl;
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif
	param_transp.build(PARAM, mf_coeft, mf_coefv);
	#ifdef M3D1D_VERBOSE_
	cout << param_transp ;
	#endif
	cout<<param_transp;
};
  
  
//! Build the tissue boundary
void
transport3d1d::build_tissue_boundary (void) 
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building tissue boundary ..." << endl;
	#endif
	BCt_transp.clear();
	BCt_transp.reserve(2*DIMT);
	// Parse BC data
	string label_in = PARAM.string_value("BClabel_transp", "Array of tissue boundary labels");
	string value_in = PARAM.string_value("BCvalue_transp", "Array of tissue boundary values");
	vector<string> labels = split(label_in, ' ');
	vector<string> values = split(value_in, ' ');
	GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
	GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
	for (unsigned f=0; f<2*DIMT; ++f) {
		BCt_transp.emplace_back(labels[f], std::stof(values[f]), 0, f);
		#ifdef M3D1D_VERBOSE_
		cout << "  face " << f << " : " << BCt_transp.back() << endl;
		#endif
	} 
	
	for (size_type bc=0; bc < BCt_transp.size(); bc++)
	cout<<BCt_transp[bc]<<endl;
	
	// Build mesht regions
	mesh_region border_faces;
	outer_faces_of_mesh(mesht, border_faces);

	for (mr_visitor i(border_faces); !i.finished(); ++i) {

		assert(i.is_face());

		// Unit outward normal : used to identify faces
		//! \todo Use getfem 5.0's function select_faces_of_normal?
		base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		if (gmm::abs(un[0] + 1.0) < 1.0E-7)      // back
			mesht.region(0).add(i.cv(), i.f());
		else if (gmm::abs(un[0] - 1.0) < 1.0E-7) // front
			mesht.region(1).add(i.cv(), i.f());
		else if (gmm::abs(un[1] + 1.0) < 1.0E-7) // left
			mesht.region(2).add(i.cv(), i.f());
		else if (gmm::abs(un[1] - 1.0) < 1.0E-7) // right
			mesht.region(3).add(i.cv(), i.f());
		else if (gmm::abs(un[2] + 1.0) < 1.0E-7) // bottom
			mesht.region(4).add(i.cv(), i.f());
		else if (gmm::abs(un[2] - 1.0) < 1.0E-7) // top
			mesht.region(5).add(i.cv(), i.f());
				
}

}

//! Build the vessel buondary
void 
transport3d1d::build_vessel_boundary(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
	#endif
try {

	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector extrema;   // global idx of extreme vertices in meshv

	Jv.clear();
	nb_extrema=0; 
	nb_junctions=0;
	
	size_type fer = nb_branches; // first empty region
	GMM_ASSERT1(meshv.has_region(fer)==0, 
		"Overload in meshv region assembling!");
	
	// List all the convexes
	dal::bit_vector nn = meshv.convex_index();
	bgeot::size_type cv;
	for (cv << nn; cv != bgeot::size_type(-1); cv << nn) {
		
		bgeot::pconvex_structure cvs = meshv.structure_of_convex(cv);
		if (cvs->nb_points()>2) 
			cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
		if (cvs->nb_faces()>2)  
			cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

		// Build regions for BCs and junctions
		// Global idx of mesh vertices
		size_type i0 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
		size_type i1 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];
		// Identify vertex type
		if (meshv.convex_to_point(i0).size()==1){ /* inflow extremum */
			// Update information
			extrema.add(i0);
			nb_extrema++;
			// Build a new region made by a single face
			GMM_ASSERT1(meshv.has_region(fer)==0, 
				"Overload in meshv region assembling!");
			meshv.region(fer).add(cv, 1);
			// Store the current index and then update it
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv_transp.size())) {
				found = (i0 == BCv_transp[bc].idx);
				if (!found) bc++;
			}
			GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");
			BCv_transp[bc].rg = fer; 
			fer++;
			// Store the containing branch index
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			BCv_transp[bc].branches.emplace_back(branch); 
		}
		else if (meshv.convex_to_point(i0).size()==2){ /* trivial inflow junction */
			// DO NOTHING
		}
		else if (meshv.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i0);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i0);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 1); // single-face region
				// Create a new junction node
				Jv.emplace_back("JUN", 0, i0, fer);
				fer++;
			}
			// Search for index of containing branch (\mathcal{P}^{in}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			// Add the inflow branch (to the right junction node)
			size_type jj = 0;
			bool found = false;
			while (!found && jj < nb_junctions){
				found = (i0 == Jv[jj].idx);
				if (!found) jj++;
			}
			//cout << "Branch -" << branch << " added to junction " << jj << endl;
			Jv[jj].value += param.R(mimv, branch);
			Jv[jj].branches.emplace_back(-branch);
			GMM_ASSERT1(branch>0, 
				"Error in network labeling: -0 makes no sense");
		}
		
		if (meshv.convex_to_point(i1).size()==1){ 
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv_transp.size())) {
				found = (i1 == BCv_transp[bc].idx);
				if (!found) bc++;
			}
			if (found){ /* outlow extremum */
				extrema.add(i1); 
				nb_extrema++; 
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				// Store the current index and then update it
				BCv_transp[bc].value *= +1.0;
				BCv_transp[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_transp[bc].branches.emplace_back(branch); 
			}
		}
		else if (meshv.convex_to_point(i1).size()==2){ /* trivial outflow junction */

			// Search for index of first containing branch (\mathcal{P}^{out}_j)
			size_type firstbranch = 0; 
			bool contained = false;
			while (!contained && firstbranch<nb_branches ) {
				contained = meshv.region(firstbranch).is_in(cv);
				if (!contained) firstbranch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i1!");

			// Check if i1 is a trivial junction (or a INT point)
			size_type cv1 = meshv.convex_to_point(i1)[0];
			size_type cv2 = meshv.convex_to_point(i1)[1];
			bool is_junc = (meshv.region(firstbranch).is_in(cv1) < 1 ||
							meshv.region(firstbranch).is_in(cv2) < 1 );
							
			if (is_junc){
				cout << "Found a trivial junction at i1 = " << i1 << endl;
				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i1);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i1);
					nb_junctions++;
					GMM_ASSERT1(meshv.has_region(fer)==0, 
						"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshv.region(fer).add(cv, 0);
					// Create a new junction node
					Jv.emplace_back("JUN", 0, i1, fer);
					fer++;
				}
				// Search for index of second containing branch (\mathcal{P}^{out}_j)
				size_type secondbranch = firstbranch+1; 
				size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
				contained = false;
				while (!contained && secondbranch<nb_branches ) {
					contained = meshv.region(secondbranch).is_in(secondcv);
					if (!contained) secondbranch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				// Add the two branches
				Jv.back().branches.emplace_back(+firstbranch);
				Jv.back().branches.emplace_back(-secondbranch);
				Jv.back().value += param.R(mimv, firstbranch);
				Jv.back().value += param.R(mimv, secondbranch);
			}
		}
		else if (meshv.convex_to_point(i1).size()>=2){ /* non-trivial outflow junction */

			// Search for index of containing branch (\mathcal{P}^{out}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");

			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i1);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i1);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 0);
				// Create a new junction node
				Jv.emplace_back("JUN", 0, i1, fer);
				// Add the outflow branch
				Jv.back().branches.emplace_back(+branch);
				Jv.back().value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << i1 << endl;
				fer++;
			}
			else {
				// Add the outflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i1 == Jv[jj].idx);
					if (!found) jj++;
				}
				Jv[jj].branches.emplace_back(+branch);
				Jv[jj].value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << jj << endl;
			}
		}

	} /* end of convexes loop */
	
	// Ckeck network assembly
	#ifdef M3D1D_VERBOSE_
	cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
	cout << "  Branches:   " << nb_branches << endl
		 << "  Vertices:   " << nn.size()+1 << endl;
	cout << "  Extrema:    " << extrema << endl;	  
	for (size_type i=0; i<BCv_transp.size(); ++i)
		cout << "    -  label=" << BCv_transp[i].label 
			 << ", value=" << BCv_transp[i].value << ", ind=" << BCv_transp[i].idx 
			 << ", rg=" << BCv_transp[i].rg << ", branches=" << BCv_transp[i].branches << endl; 
	cout << "  Junctions: " << junctions << endl;
	for (size_type i=0; i<Jv.size(); ++i)
		cout << "    -  label=" << Jv[i].label 
			 << ", value=" << Jv[i].value << ", ind=" << Jv[i].idx 
			 << ", rg=" << Jv[i].rg << ", branches=" << Jv[i].branches << endl; 
	cout << "---------------------------------------- "   << endl;
	#endif

} 
GMM_STANDARD_CATCH_ERROR; // catches standard errors

} /* end of build_vessel_boundary */


//! Assembly the transport problem
void transport3d1d::assembly (void)
{
 	std::cout<<"assemble transport problem"<<std::endl;
 	//1 Build the monolithic matrix AM
	assembly_mat(); 
 	} // end of assembly
 
//! Assembly the transport matrices
void 
transport3d1d::assembly_mat(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM, UM, FM ..." << endl;
	#endif
	gmm::resize(AM_transp_t, dof_transp.Ct(), dof_transp.Ct()); gmm::clear(AM_transp_t);
	gmm::resize(UM_transp_t, dof_transp.Ct()); gmm::clear(UM_transp_t);
	gmm::resize(FM_transp_t, dof_transp.Ct()); gmm::clear(FM_transp_t);
	gmm::resize(AM_transp_v, dof_transp.Cv(), dof_transp.Cv()); gmm::clear(AM_transp_v);
	gmm::resize(UM_transp_v, dof_transp.Cv()); gmm::clear(UM_transp_v);
	gmm::resize(FM_transp_v, dof_transp.Cv()); gmm::clear(FM_transp_v);
	#ifdef M3D1D_VERBOSE_
	cout << "the monolithic matrix AM ..." << endl;
	#endif
	// Reaction matrix for the interstitial problem
	sparse_matrix_type Rt(dof_transp.Ct(), dof_transp.Ct()); gmm::clear(Rt);
	// Diffusion matrix for the interstitial problem
	sparse_matrix_type Dt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Dt);
	//Transport matrix for interstitial problem
	sparse_matrix_type Bt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Bt);
	// Mass(time derivative)  matrix for the interstitial problem
	sparse_matrix_type Mt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Mt);
	// Adhesion  matrix for the network problem					
	sparse_matrix_type Adhv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Adhv);	
		
	// Reaction matrix for the network problem
	sparse_matrix_type Rv(dof_transp.Cv(), dof_transp.Cv()); gmm::clear(Rv);
	// Diffusion matrix for the network problem
	sparse_matrix_type Dv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Dv);
	//Transport matrix for network problem
	sparse_matrix_type Bv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bv);
	// Mass (time derivative)  matrix for the network problem
	sparse_matrix_type Mv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mv);

	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Btt);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof_transp.Ct(), dof_transp.Cv());gmm::clear(Btv);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bvv);
	
	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mbar);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mlin);
	
	sparse_matrix_type Btt_(dof.Pt(), dof.Pt()); gmm::clear(Btt_);
	sparse_matrix_type Btv_(dof.Pt(), dof.Pv()); gmm::clear(Btv_);
	sparse_matrix_type Bvt_(dof.Pv(), dof.Pt()); gmm::clear(Bvt_);
	sparse_matrix_type Bvv_(dof.Pv(), dof.Pv()); gmm::clear(Bvv_);
	sparse_matrix_type Mbar_(dof.Pv(), dof.Pt()); gmm::clear(Mbar_);
	sparse_matrix_type Mlin_(dof.Pv(), dof.Pt()); gmm::clear(Mlin_);

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Rt, Mt and Dt ..." << endl;
	#endif
	
	// Build the reaction coefficient:
	// The reaction Dalpha is not considered (NPs do not degrade);
	// The lymphatic term (Q_pl(Pt-Pl)ct,bt)_\Omega does not appear -> it cancels out with the non conservative term of advection (\nabla \cdot u_t ct, bt )
	// indeed from the flow equation (divergence) : \nabla \cdot u_t = - Q_pl(Pt-Pl) +(...)\delta_Lambda -> neglecting the delta term

	vector_type mass_coeff(dof.Pt()); gmm::clear(mass_coeff);
	vector_type Pl(dof.Pt(),PARAM.real_value("PL")); 
	gmm::scale(Pl, -1.0); 
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut(), dof.Pt())) ,  mass_coeff);
	gmm::add(Pl ,  mass_coeff);
	gmm::scale (param_transp.Q_pl(), mass_coeff);
	gmm::add(param_transp.Dalpha(), mass_coeff);  



	/* Tissue artificial diffusion (if advection)
	modify the diffusion coefficient:
	At_stab=At(1+Pe_loc); with Pe_loc=Ut_dim * h/ 2*Dt= Ut_adim * U * h * d / 2 * At * U * d = Ut_adim * h / 2 * At;
	*/
	vector_type element_size_t(mesht.dim());
	size_type kk=0;
	for(size_type i=0; i < mesht.dim(); i++){
		element_size_t[kk]=mesht.convex_area_estimate(i,2);
		kk++;
		}
	scalar_type max_size_t=*max_element(element_size_t.begin(), element_size_t.end());

	//find the maximum of the velocity field
	scalar_type max_U_t=0.0;
	scalar_type max_U_positive_t=0.0;
	scalar_type max_U_negative_t=0.0;
	
	vector_type Ut_diff_art(dof.Ut()); gmm::clear(Ut_diff_art);
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())) ,  Ut_diff_art);

	max_U_positive_t=*max_element(Ut_diff_art.begin(), Ut_diff_art.end());
	max_U_negative_t=*min_element(Ut_diff_art.begin(), Ut_diff_art.end());

	if(max_U_positive_t > fabs(max_U_negative_t))
		max_U_t=max_U_positive_t;
	else
		max_U_t=fabs(max_U_negative_t);
	

	scalar_type artif_diff_t=param_transp.At(1)+(max_U_t*max_size_t/2.0);

	vector_type At_stabilized(mf_coeft.nb_dof(),artif_diff_t);

	bool ADVECTION_T = PARAM.int_value("ADVECTION_T", "flag for advection term in tissue");
	if(ADVECTION_T ==1){ // Artificial diffusion
		asm_tissue_darcy_transp(Rt, Dt, Mt, mimt, mf_Ct, mf_Pt, mf_coeft, mass_coeff, At_stabilized ); 
	}
	if(ADVECTION_T ==0){ //NO artificial diffusion
		asm_tissue_darcy_transp(Rt, Dt, Mt, mimt, mf_Ct, mf_Pt, mf_coeft, mass_coeff, param_transp.At() ); 
	}

	gmm::scale(Mt, (1.0/param_transp.dt())); // dt time step

	
	
	// Copy Rt: reaction term
	bool REACTION = PARAM.int_value("REACTION", "flag for reaction term");
	if(REACTION ==1){
	gmm::add(Rt, AM_transp_t); 
	}
	 
	// Copy Mt: time 
	if(descr_transp.STATIONARY ==0)
	gmm::add(Mt, AM_transp_t); 


	// Copy Dt: diffusion
	bool DIFFUSION_T = PARAM.int_value("DIFFUSION_T", "flag for diffusion term in tissue");
	if(DIFFUSION_T ==1){	
	gmm::add(Dt, AM_transp_t); 
	}	

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mv and Dv ..." << endl;
	#endif




	/* Network artificial diffusion (if advection)
	modify the diffusion coefficient:
	Av_stab=Av(1+Pe); con Pe_loc=Uv_dim * h/ 2*Dv= Uv_adim * U * h * d / 2 * Av * U * d = Uv_adim * h / 2 * Av;
	*/

	vector_type element_size(meshv.dim());
	size_type j=0;
	for(size_type i=0; i < meshv.dim(); i++){
		element_size[j]=meshv.convex_area_estimate(i,2);
		j++;
		}
	scalar_type max_size=*max_element(element_size.begin(), element_size.end());

	//find the maximum of the velocity field
	scalar_type max_U=0.0;
	size_type shift3 =0;	
	vector_type max_U_pos_vec(nb_branches);gmm::clear(max_U_pos_vec);
	vector_type max_U_neg_vec(nb_branches);gmm::clear(max_U_neg_vec);	
	scalar_type max_U_positive=0.0;
	scalar_type max_U_negative=0.0;
	
	for (size_type i=0; i<nb_branches; ++i){
		if(i>0) shift3 += mf_Uvi[i-1].nb_dof();
		vector_type Uvi2(mf_Uvi[i].nb_dof()); gmm::clear(Uvi2);
		gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift3, mf_Uvi[i].nb_dof())) ,  Uvi2);
		max_U_pos_vec[i]=*max_element(Uvi2.begin(), Uvi2.end()); //max on a branch
		max_U_neg_vec[i]=*min_element(Uvi2.begin(), Uvi2.end()); //min on a branch
 	}
	
	max_U_positive=*max_element(max_U_pos_vec.begin(), max_U_pos_vec.end());
	max_U_negative=*min_element(max_U_neg_vec.begin(), max_U_neg_vec.end());

	if(max_U_positive > fabs(max_U_negative))
		max_U=max_U_positive;
	else
		max_U=fabs(max_U_negative);
	
	
	scalar_type artif_diff=param_transp.Av(1)+(max_U*max_size/2.0);
	vector_type Av_stabilized(mf_coefv.nb_dof(),artif_diff);

	// Build Mvvi and Dvvi
	bool ADVECTION_V = PARAM.int_value("ADVECTION_V", "flag for advection term in vessels");
	if(ADVECTION_V ==1){  // Artificial diffusion
		asm_network_poiseuille_transp(Dv, Mv, mimv,mf_Cv, mf_coefv, Av_stabilized);
	}
	if(ADVECTION_V ==0){  // NO Artificial diffusion
		asm_network_poiseuille_transp(Dv, Mv, mimv,mf_Cv, mf_coefv, param_transp.Av());
	}
	gmm::scale(Mv, (1.0/param_transp.dt()));
				
	// Copy Mv: time 
	if(descr_transp.STATIONARY ==0)
	gmm::add(Mv, AM_transp_v); 
		
	// Copy Dv: diffusion 
	bool DIFFUSION_V = PARAM.int_value("DIFFUSION_V", "flag for diffusion term in vessel");
	if(DIFFUSION_V ==1){	
	gmm::add(Dv, AM_transp_v);
	}			

		

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling the adhesive term ..." << endl;
	#endif	
	
	// Build the adhesion term 
	
	//Create the vectors for the projection on the coefv_mesh (global on the network)

	gmm::resize(Reproj, mf_coefv.nb_dof()); gmm::clear(Reproj);
	gmm::resize(WSRproj, mf_coefv.nb_dof()); gmm::clear(WSRproj);
	gmm::resize(P_adhproj, mf_coefv.nb_dof()); gmm::clear(P_adhproj);
	gmm::resize(PiGrecoproj, mf_coefv.nb_dof()); gmm::clear(PiGrecoproj);
	gmm::resize(PiGrecoproj_scaled, mf_coefv.nb_dof()); gmm::clear(PiGrecoproj_scaled);	
	//We comment them becacause we don't project PiGreco on Cv anymore
	//gmm::resize(PiGrecoproj_cv,dof_transp.Cv()); gmm::clear(PiGrecoproj_cv);
	//gmm::resize(PiGrecoproj_scaled_cv, dof_transp.Cv()); gmm::clear(PiGrecoproj_scaled_cv);

	//mass matrix for the coefficient
	gmm::resize (Mcc,mf_coefv.nb_dof(), mf_coefv.nb_dof()); gmm::clear(Mcc);
	getfem::asm_mass_matrix(Mcc, mimv, mf_coefv);
	//mass matrix for the coefficient
	//sparse_matrix_type Mcc_cv(dof_transp.Cv(), dof_transp.Cv()); gmm::clear(Mcc_cv);
	//getfem::asm_mass_matrix(Mcc_cv, mimv, mf_Cv);

	//rhs term
	vector_type f_Re(mf_coefv.nb_dof()); gmm::clear(f_Re);
	vector_type f_WSR(mf_coefv.nb_dof()); gmm::clear(f_WSR);
	vector_type f_Padh(mf_coefv.nb_dof()); gmm::clear(f_Padh);
	vector_type f_Pigreco(mf_coefv.nb_dof()); gmm::clear(f_Pigreco);
	vector_type f_Pigreco_scaled(mf_coefv.nb_dof()); gmm::clear(f_Pigreco_scaled);
	//vector_type f_Pigreco_cv(dof_transp.Cv()); gmm::clear(f_Pigreco_cv);
	//vector_type f_Pigreco_scaled_cv(dof_transp.Cv()); gmm::clear(f_Pigreco_scaled_cv);
	
	bool NANO = PARAM.int_value("NANO", "flag for adhesive term in vessel");
	bool SATURATION  = PARAM.int_value("SATURATION","Flag to use the saturation per the adhesive term");

	size_type shift1 = 0; //dof of mf_Uvi[i]
	size_type shift2 = 0; //dof of mf_coefvi[i]
	for (size_type i=0; i<nb_branches; ++i){

		// Wall shear rate vector on a branch
		vector_type WSRi (mf_Uvi[i].nb_dof()); gmm::clear(WSRi); 
		// Probability of adhesion vector on a branch
		vector_type P_adhi (mf_Uvi[i].nb_dof()); gmm::clear(P_adhi); 
		// Reynolds vector on a branch
		vector_type Rei (mf_Uvi[i].nb_dof()); gmm::clear(Rei);
	
		if(i>0) shift1 += mf_Uvi[i-1].nb_dof();
		if(i>0) shift2 += mf_coefvi[i-1].nb_dof();
		// Velocity of the i-th branch
 		vector_type Uvi1(mf_Uvi[i].nb_dof()); gmm::clear(Uvi1);
		vector_type Uvi1_abs(mf_Uvi[i].nb_dof()); gmm::clear(Uvi1_abs);

		gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift1, mf_Uvi[i].nb_dof())) ,  Uvi1);


		for (size_type k=0; k<mf_Uvi[i].nb_dof(); ++k){
			Uvi1_abs[k]=fabs(Uvi1[k]);
		}
		

		// Radius of the i-th branch
		scalar_type Ri = compute_radius(mimv, mf_coefv, param.R(), i); ; 
		
		//Wall shear rate		
		wall_shear_stress(WSRi, Ri, Uvi1_abs); //Adimensional for the export 
		vector_type WSRi_dim(mf_Uvi[i].nb_dof()); gmm::clear(WSRi_dim); //to be put in P_a
		gmm::add(WSRi,WSRi_dim);		
		gmm::scale(WSRi_dim, (param_transp.Uadim()/ param_transp.dadim()));// Rescale WSRi_adim * U / d -> WSR is dimensional

		// Reynolds
		reynolds (Rei, Uvi1_abs, Ri, param.mu_v(), param_transp.rho()); // Dimensional radius and velocity 
		gmm::scale(Rei,(param_transp.Uadim()*param_transp.dadim())); //Rescale Re_adimensi * U * d -> Rei is the real Reynolds
		
		//Probability of adhesion
		bool PA_INPUT  = PARAM.int_value("PA_INPUT");
		//If probability of adhesion from the exponential function
		if(PA_INPUT==0){
		scalar_type alpha2=param_transp.m_r()*(1.0-pow(1.0-(2.0*param_transp.h_0()/param_transp.dp()),2.0));
		probability_adhesion0(P_adhi, WSRi_dim, param_transp.m_l(), param_transp.Ka(),alpha2,param_transp.r_0(),param_transp.beta_nano(),param.mu_v(),param_transp.m_r()); //Dimensional WRS
		}
		// If probability of adhesion from Lattice Boltzmann method
		// Case strong bond -> if possible: modify this and pass these values as input
		if(PA_INPUT==1){
		scalar_type Pa_min;
		scalar_type Pa_med;
		scalar_type Pa_max;
			if(param_transp.rho_l()==0.3){Pa_min=0.3; Pa_med=0.0; Pa_max=0.0;}
			else if (param_transp.rho_l()==0.5){Pa_min=1.0; Pa_med=0.33; Pa_max=0.28;}
			else if (param_transp.rho_l()==0.7){Pa_min=1.0; Pa_med=0.64; Pa_max=0.56;}
			else if (param_transp.rho_l()==0.9){Pa_min=1.0; Pa_med=0.8; Pa_max=0.7;}
			else {GMM_ASSERT1(PA_INPUT==1, "rho_l value not valid");}
		probability_adhesion(P_adhi, Rei, Pa_min, Pa_med, Pa_max);	
		}
		
		
	
		vector_type adhesion_coeff(mf_Uvi[i].nb_dof()); gmm::clear(adhesion_coeff); 
		vector_type adhesion_coeff_scaled(mf_Uvi[i].nb_dof()); gmm::clear(adhesion_coeff_scaled); 
		//\Pi(s) = P_adh* WSR *dp/2
		vector_type temporary(mf_Uvi[i].nb_dof()); gmm::clear(temporary); 
		gmm::add(WSRi_dim, temporary); 			// P_adh * WSR (component by component) -> P_adh; 
		gmm::scale(P_adhi, temporary); 			// P_adh*WSR -> adhesion
		gmm::scale(temporary, (param_transp.dp()/2.0)); // Scale dp/2
		gmm::scale(temporary, (1.0/param_transp.Uadim())); // Adimensionalization
		gmm::add(temporary,adhesion_coeff);
		gmm::add(temporary,adhesion_coeff_scaled);
		gmm::scale(adhesion_coeff_scaled, 2.0/Ri); 
		
		if(NANO==1 && SATURATION==0){
		asm_network_nano_transp(Adhv, mimv, mf_Cv, mf_Uvi[i], adhesion_coeff_scaled, meshv.region(i));	
		masslumping(Adhv); 
		}

		//Build the rhs terms for the projection
		getfem::generic_assembly assem;
		assem.push_mi(mimv);
		assem.push_mf(mf_Uvi[i]);
		assem.push_mf(mf_coefv);
		//assem.push_mf(mf_Cv);
		assem.push_data(Rei);
		assem.push_data(WSRi);
		assem.push_data(P_adhi);
		assem.push_data(adhesion_coeff);
		assem.push_data(adhesion_coeff_scaled);
		assem.push_vec(f_Re);
		assem.push_vec(f_WSR);
		assem.push_vec(f_Padh);
		assem.push_vec(f_Pigreco); 
		assem.push_vec(f_Pigreco_scaled);		// Project Pi on coefv(PK(1,0)) for export_vtk
		//assem.push_vec(f_Pigreco_cv); 		// Project Pi on cv(PK(1,1)) for computing Psi
		//assem.push_vec(f_Pigreco_scaled_cv); 	// Project Pi*2/R on cv(PK(1,1)) for Pigrecostar
	assem.set("m1=data$1(#1);"
			  "V$1(#2)+=comp(Base(#1).Base(#2))(j,:).m1(j);"
			  "m2=data$2(#1);"
			  "V$2(#2)+=comp(Base(#1).Base(#2))(j,:).m2(j);"
			  "m3=data$3(#1);"
			  "V$3(#2)+=comp(Base(#1).Base(#2))(j,:).m3(j);"
			  "m4=data$4(#1);"
			  "V$4(#2)+=comp(Base(#1).Base(#2))(j,:).m4(j);"
			  "m5=data$5(#1);"
			  "V$5(#2)+=comp(Base(#1).Base(#2))(j,:).m5(j);"
			 // "m6=data$4(#1);"
			 // "V$6(#3)+=comp(Base(#1).Base(#3))(j,:).m6(j);"
			); 
		assem.assembly(meshv.region(i));

		
	}
	
	// Solve the projection problem
	#ifdef M3D1D_VERBOSE_
	cout << "  Solve the projection of the fluid dynamical quantities ... " << endl;
	#endif
	gmm::csc_matrix<scalar_type> Mcc_aux;
	gmm::clean(Mcc, 1E-12);
	gmm::copy(Mcc, Mcc_aux);
	/*gmm::csc_matrix<scalar_type> Mcc_cv_aux;
	gmm::clean(Mcc_cv, 1E-12);
	gmm::copy(Mcc_cv, Mcc_cv_aux);
	*/
	scalar_type cond;

	//Reynolds
	gmm::SuperLU_solve(Mcc_aux, Reproj, f_Re,cond);

	//WSR
	gmm::SuperLU_solve(Mcc_aux, WSRproj, f_WSR,cond);

	//P_adhesion
	gmm::SuperLU_solve(Mcc_aux, P_adhproj, f_Padh,cond);
		
	//Vascular adhesion parameter PiGreco (adimensional) on coefv--> PiGreco/Uadim
	gmm::SuperLU_solve(Mcc_aux, PiGrecoproj, f_Pigreco,cond);
		
	//Vascular adhesion parameter PiGreco*2/R (adimensional) on coefv--> 2/R*PiGreco/Uadim
	gmm::SuperLU_solve(Mcc_aux, PiGrecoproj_scaled, f_Pigreco_scaled,cond);
	
	//Vascular adhesion parameter PiGreco (adimensional) on base cv--> PiGreco/Uadim
	//gmm::SuperLU_solve(Mcc_cv_aux, PiGrecoproj_cv, f_Pigreco_cv,cond);
	
	//Vascular adhesion parameter PiGreco (adimensional) on cv--> PiGreco/Uadim
	//gmm::SuperLU_solve(Mcc_cv_aux, PiGrecoproj_scaled_cv, f_Pigreco_scaled_cv,cond);
	

	//gmm::clear(Mcc); //Don't clear Mcc beacuse we need it after for projecting Cvold
	//gmm::clear(Mcc_cv);
	gmm::clear(f_Re);
	gmm::clear(f_WSR);
	gmm::clear(f_Padh);
	gmm::clear(f_Pigreco);
	gmm::clear(f_Pigreco_scaled);
	//gmm::clear(f_Pigreco_cv);
	//gmm::clear(f_Pigreco_scaled_cv);


	if(NANO==1 && SATURATION==0){
	gmm::add(Adhv, AM_transp_v);
	}
		

	//ADVECTION	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Bt and Bv ..." << endl;
	#endif		

	if(ADVECTION_T ==1){
	// advection tissue term: 				
	vector_type Ut(dof.Ut()); gmm::clear(Ut);
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())) ,  Ut);
	asm_advection_tissue(Bt, mimt, mf_Ct, mf_Ut,Ut);

	gmm::add(Bt, AM_transp_t); 		
		
	}				

  
	
	if(ADVECTION_V ==1){
	// advection vessel term: 
	size_type shift = 0; // initial position of each branch
	
	//Tangent vectors lambda
	vector_type lambdax; // tangent versor: x component
	vector_type lambday; // tangent versor: y component
	vector_type lambdaz; // tangent versor: z component
	std::ifstream ifs(descr.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr.MESH_FILEV);
	asm_tangent_versor(ifs, lambdax, lambday, lambdaz);
	ifs.close();
	
	
	//loop on each branch
	for (size_type i=0; i<nb_branches; ++i){

	string nramo = "";
	std::ostringstream convert;
	convert << i;
	nramo = convert.str();
	
	// Initial position of each branch
	if(i>0) shift += mf_Uvi[i-1].nb_dof();
	// Velocity Uvi of the i-th branch
 	vector_type Uvi(mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);
	vector_type Uvi_abs(mf_Uvi[i].nb_dof()); gmm::clear(Uvi_abs);
	for(size_type k=0; k<mf_Uvi[i].nb_dof(); k++){
		Uvi_abs[k]=fabs(Uvi[k]);
	}
	// Tangent versors of the i-th branch
	vector_type lambdax_K, lambday_K, lambdaz_K;
	for(size_type j=0; j<mf_coefvi[i].nb_dof(); j++){
			lambdax_K.emplace_back(lambdax[i]);
			lambday_K.emplace_back(lambday[i]);
			lambdaz_K.emplace_back(lambdaz[i]);
	}

	// Advection matrix
	asm_advection_network(Bv, mimv, mf_Cv, mf_coefvi[i], mf_Uvi[i], Uvi, lambdax_K, lambday_K, lambdaz_K, meshv.region(i) );


	} //end of loop on the branches

	gmm::add(Bv, AM_transp_v);
	

	//Build Rv: assemble int(c_v* d(u_v)/ds, phi_i) = int(c_v* [Bvv*Pv-Bvt*Pt], phi_i) --> for the divergence of uv
	

	asm_exchange_aux_mat(Mbar_, Mlin_, mimv, mf_Pt, mf_Pv, param.R(), descr.NInt);
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION");
	asm_exchange_mat(Btt_, Btv_, Bvt_, Bvv_, mimv, mf_Pv, mf_coefv, Mbar_, Mlin_, param.Q(), NEWFORM);

	vector_type Rv_coef (dof.Pv());
	gmm::mult(Bvv_, 
		  gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
		  Rv_coef);
		  	  
	gmm::mult_add(gmm::scaled(Bvt_, -1.0),
		      gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut(), dof.Pt())),
		      Rv_coef);		
	// Oncotic term
	scalar_type Pi_t=PARAM.real_value("Pi_t", "Interstitial Oncotic Pressure [-]");
	scalar_type Pi_v=PARAM.real_value("Pi_v", "Plasmatic Oncotic Pressure [-]");
	scalar_type sigma=PARAM.real_value("sigma", "Reflection Coefficient sigma");
	scalar_type picoef=sigma*(Pi_v-Pi_t);
        vector_type DeltaPi(dof.Pv(),picoef);
        gmm::mult_add(scaled(Bvv_,-1.0),DeltaPi,Rv_coef);
			
	// Build Rv				
	getfem::asm_mass_matrix_param(Rv, mimv, mf_Cv, mf_Pv, Rv_coef);
	// Copy Rv
	gmm::add(gmm::scaled(Rv,-1.0/(pi*param.R(0)*param.R(0))), AM_transp_v);

	}


	// De-allocate memory
	gmm::clear(Mt);  gmm::clear(Dt); 
	gmm::clear(Mv); gmm::clear(Dv);
	gmm::clear(Bt);  gmm::clear(Bv);
	gmm::clear(Mbar);  gmm::clear(Mlin);
	gmm::clear(Btt);  gmm::clear(Btv);
	gmm::clear(Bvt);  gmm::clear(Bvv);
	gmm::clear(Btt_);  gmm::clear(Btv_);
	gmm::clear(Bvt_);  gmm::clear(Bvv_);
	gmm::clear(Mbar_);  gmm::clear(Mlin_);
	gmm::clear(Adhv);

}

void 
transport3d1d::assembly_rhs_tissue(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term in the RHS..." << endl;
	#endif

	sparse_matrix_type Att(dof_transp.Ct(),   dof_transp.Ct());
	vector_type Ft(dof_transp.Ct());
	
	gmm::add(AM_temp_t, Att);
	gmm::scale(AM_temp_t, 0.0);	
				
	gmm::add(FM_temp_t,Ft);	 
	gmm::scale(FM_temp_t, 0.0);

	
	scalar_type beta_t = PARAM.real_value("BETAtissue_transp", "Coefficient for mixed BC for transport problem in tissue");
	asm_tissue_bc_transp(Att, Ft, mimt, mf_Ct, mf_coeft, BCt_transp, beta_t);
	gmm::add(Att, AM_temp_t);
	gmm::add(Ft, FM_temp_t);

	// De-allocate memory
	gmm::clear(Att);
	gmm::clear(Ft);	
}

void 
transport3d1d::assembly_rhs_network(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term  in the RHS..." << endl;
	#endif
		
	// Build the boundary conditions: for AM_temp and Fv

	sparse_matrix_type Avv(dof_transp.Cv(), dof_transp.Cv());
	vector_type Fv(dof_transp.Cv());
	gmm::add(AM_temp_v, Avv);
	gmm::scale(AM_temp_v, 0.0);	
				
	gmm::add(FM_temp_v, Fv);	
	gmm::scale(FM_temp_v, 0.0);

	scalar_type beta_v  = PARAM.real_value("BETAvessel_transp", "Coefficient for mixed BC for transport problem in vessels");
		cout<< "inj_time....."<<param_transp.inj_time()<<endl;
		if(t<=param_transp.inj_time())
	asm_network_bc_transp(Avv, Fv, mimv, mf_Cv, mf_coefv, BCv_transp, beta_v );
//===========================================================================	
	else{
	//vector< node > BCv_transp_noinjection; gmm::clear( BCv_transp_noinjection);
	for(size_type i=0; i<BCv_transp.size(); i++){
		//BCv_transp_noinjection[i].label = BCv_transp[i].label;
		//BCv_transp_noinjection[i].rg = BCv_transp[i].rg;
		if(BCv_transp[i].label=="DIR")
			//BCv_transp_noinjection[i].value = 0;
			BCv_transp[i].value = 0;
		//else
		//	BCv_transp_noinjection[i].value = BCv_transp[i].value;
	}
	
	//asm_network_bc_transp(Avv, Fv, mimv, mf_Cv, mf_coefv, BCv_transp_noinjection, beta_v );
	asm_network_bc_transp(Avv, Fv, mimv, mf_Cv, mf_coefv, BCv_transp, beta_v );
	}
//===========================================================================
	gmm::add(Avv, AM_temp_v);
	gmm::add(Fv, FM_temp_v);

	// De-allocate memory
	gmm::clear(Avv);
	gmm::clear(Fv);

	
}

void transport3d1d::update_tissue (vector_type RHS_nano){ 

	// Assembled AM matrix (Do not modify the time independent terms - Modify the time dep terms)
	// Empty FM_temp rhs (Add time dependence; 3D BCs)
	gmm::copy(AM_transp_t, AM_temp_t);
	gmm::copy(FM_transp_t, FM_temp_t);


	// Temporal update: F + F_old

	// update rhs (there is the time step mass term)
	vector_type TFt(dof_transp.Ct());
	asm_source_term(TFt, mimt, mf_Ct, mf_Ct, UM_transp_t);
	gmm::scale(TFt, (1.0/param_transp.dt())); // dt time step
	gmm::add(TFt, FM_temp_t);
	
	// add to the rhs the source term
	gmm::scale(RHS_nano, 2*pi*param.R(0)/(param_transp.cNP()*param_transp.Uadim())); 	
	gmm::add(RHS_nano, FM_temp_t);	

	gmm::clear(UM_transp_t);
	gmm::clear(TFt);

	assembly_rhs_tissue();


}

void transport3d1d::update_network (vector_type Pigreco){

	// Assembled AM matrix (Do not modify the time independent terms - Modify the time dep terms)
	// Empty FM_temp rhs (Add time dependence; 3D BCs)
	gmm::copy(AM_transp_v, AM_temp_v);
	gmm::copy(FM_transp_v, FM_temp_v);

	bool SATURATION  = PARAM.int_value("SATURATION","Flag to use the saturation per the adhesive term");

	if(SATURATION==1){
		sparse_matrix_type Adhv(dof_transp.Cv(), dof_transp.Cv()); gmm::clear(Adhv);
		asm_network_nano_transp(Adhv, mimv, mf_Cv, mf_coefv, Pigreco);	
		masslumping(Adhv);
		gmm::add(Adhv, AM_temp_v);
	}



	// Temporal update: F + F_old
	// update rhs (there is the time step mass term)

	vector_type TFv(dof_transp.Cv());
	asm_source_term(TFv,mimv, mf_Cv, mf_Cv, UM_transp_v);
	gmm::scale(TFv, (1.0/param_transp.dt())); // dt time step

	gmm::add(TFv, FM_temp_v);
	gmm::clear(UM_transp_v);
	gmm::clear(TFv);

	assembly_rhs_network();


}

 bool transport3d1d::solve (void)
 {
	#ifdef M3D1D_VERBOSE_
  	std::cout<<"Solving transport problem"<<std::endl;
	#endif

	gmm::resize(AM_temp_t, dof_transp.Ct(), dof_transp.Ct()); gmm::clear(AM_temp_t);
	gmm::resize(FM_temp_t, dof_transp.Ct()); gmm::clear(FM_temp_t);
	gmm::resize(AM_temp_v, dof_transp.Cv(), dof_transp.Cv()); gmm::clear(AM_temp_v);
	gmm::resize(FM_temp_v, dof_transp.Cv()); gmm::clear(FM_temp_v);
	
	double time = gmm::uclock_sec();
	double time_count = 0;	
	int iteraz = 0;
	
	bool SATURATION  = PARAM.int_value("SATURATION","Flag to use the saturation per the adhesive term");
	bool NANO  = PARAM.int_value("NANO");

	#ifdef M3D1D_VERBOSE_
  	std::cout<<"Solving vessel problem"<<std::endl;
	#endif

	vector_type Cvold(dof_transp.Cv()); gmm::clear(Cvold);
	vector_type Cvnew(dof_transp.Cv()); gmm::clear(Cvnew);

	vector_type Psiold(mf_coefv.nb_dof()); gmm::clear(Psiold);
	gmm::resize(Psi, mf_coefv.nb_dof()); gmm::clear(Psi);
	
	vector_type psi_max_vec(mf_coefv.nb_dof()); gmm::clear(psi_max_vec);
	vector_type Pigrecostar(mf_coefv.nb_dof()); gmm::clear(Pigrecostar);
	 
	if(NANO==1 && SATURATION==1){
	psi_max_vec.assign(mf_coefv.nb_dof(),  param_transp.psi_max());
	}



	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar1(dof_transp.Cv(), dof_transp.Ct()); gmm::clear(Mbar1);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin1(dof_transp.Cv(), dof_transp.Ct()); gmm::clear(Mlin1);
	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt1(dof_transp.Ct(), dof_transp.Ct()); gmm::clear(Btt1);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv1(dof_transp.Ct(), dof_transp.Cv()); gmm::clear(Btv1);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt1(dof_transp.Cv(), dof_transp.Ct()); gmm::clear(Bvt1);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv1(dof_transp.Cv(), dof_transp.Cv()); gmm::clear(Bvv1);

	asm_exchange_aux_mat(Mbar1, Mlin1, mimv, mf_Ct, mf_Cv, param.R(), descr.NInt);
	const bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
	
	asm_exchange_mat_transp(Btt1, Btv1, Bvt1, Bvv1, mimv, mf_Cv, Mbar1, Mlin1, NEWFORM);


	
	std::ofstream outMeanCt("MeanCt.txt"); // Average of C_t on the 3D domain at each time step

	for(t=0;t<=param_transp.T()*(!descr_transp.STATIONARY) ; t = t + param_transp.dt() + (param_transp.dt()==0) ){ 
	time_count++; 
	iteraz++; 
	std::cout<<"iteration number:"<<time_count<<std::endl;
	std::cout<<"time = "<<t<<" s"<<std::endl;	
	
	//if(NANO==1 && SATURATION==1 && t==0)
	if(NANO==1 && SATURATION==1)
	{


		vector_type f_Cvold_proj(mf_coefv.nb_dof()); gmm::clear(f_Cvold_proj);
		vector_type Cvold_proj(mf_coefv.nb_dof()); gmm::clear(Cvold_proj);
		
	    //Projection of Cvold on PK(1,0) for computing Psi	
		//Build the rhs terms for the projection
		getfem::generic_assembly assem;
		assem.push_mi(mimv);
		assem.push_mf(mf_Cv);
		assem.push_mf(mf_coefv);
		assem.push_data(Cvold);
		assem.push_vec(f_Cvold_proj); 		
		assem.set("m1=data$1(#1);"
			  "V$1(#2)+=comp(Base(#1).Base(#2))(j,:).m1(j);");
		assem.assembly();

		gmm::csc_matrix<scalar_type> Mcc_aux1;
		gmm::clean(Mcc, 1E-12);
		gmm::copy(Mcc, Mcc_aux1);
		scalar_type cond;
			
		//solve the projection
		gmm::SuperLU_solve(Mcc_aux1, Cvold_proj, f_Cvold_proj,cond);

		for(size_type i=0; i<mf_coefv.nb_dof(); i++) { 

				Psi[i] = PiGrecoproj[i]*std::max(0.0,(psi_max_vec[i]-Psiold[i]))*Cvold_proj[i];
		}
		gmm::scale(Psi, param_transp.dadim()*param_transp.dt()/param_transp.psi_max());
		gmm::add(Psiold, Psi);
		
		for(size_type i=0; i<mf_coefv.nb_dof(); i++) { 
      			Pigrecostar[i] = PiGrecoproj_scaled[i]*std::max(0.0,(psi_max_vec[i]-Psi[i]));
      		}
		gmm::scale(Pigrecostar, 1.0/param_transp.psi_max());
		
		gmm::clear(Psiold);
		gmm::add(Psi, Psiold);

				//Print Psi at several time steps
		string time_suff = "";
		std::ostringstream convert;
		convert << time_count;
		time_suff = convert.str();
		
		//if(iteraz%5==0){ 
		
		cout << "  Exporting Psi ..." << endl;
		vtk_export exp_Psi(descr_transp.OUTPUT+"Psi_t"+time_suff+".vtk");
		exp_Psi.exporting(mf_coefv);
		exp_Psi.write_mesh();
		exp_Psi.write_point_data(mf_coefv, Psi, "Psi");


		//}
		cout << "End exporting_Psi"<< endl;
	}	
	

	update_network(Pigrecostar);

	gmm::csc_matrix<scalar_type> A_transp_v;
	gmm::clean(AM_transp_v, 1E-12);
	gmm::copy(AM_temp_v, A_transp_v);
	
	vector_type F_transp_v(gmm::vect_size(FM_transp_v));
	gmm::clean(FM_transp_v, 1E-12);
	gmm::copy(FM_temp_v, F_transp_v);
	

	if ( descr_transp.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method to network problem... " << endl;
		#endif
		scalar_type cond_v;
		gmm::SuperLU_solve(A_transp_v, UM_transp_v, F_transp_v, cond_v);
		cout << "  Condition number (transport network problem): " << cond_v << endl;
	}
	/*else { // Iterative solver //

		// Iterations
		gmm::iteration iter(descr_transp.RES);  // iteration object with the max residu
		iter.set_noisy(1);               // output of iterations (2: sub-iteration)
		iter.set_maxiter(descr_transp.MAXITER); // maximum number of iterations

		// Preconditioners
		//! \todo Add preconditioner choice to param file
		// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
		gmm::identity_matrix PM; // no precond
		//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
		//gmm::ilu_precond<sparse_matrix_type> PM(AM);
		// ...
		//gmm::clear(AM);
		// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15
	
		if ( descr_transp.SOLVE_METHOD == "CG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Conjugate Gradient method ... " << endl;
			#endif
			gmm::identity_matrix PS;  // optional scalar product
			gmm::cg(AM_transp, UM_transp, F_transp, PS, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "BiCGstab" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
			#endif
			gmm::bicgstab(AM, UM, FM, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "GMRES" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Generalized Minimum Residual method ... " << endl;
			#endif
			size_type restart = 50;
			gmm::gmres(A_transp, UM, FM, PM, restart, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "QMR" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
			#endif
			gmm::qmr(AM, UM, FM, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "LSCG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
			#endif
			gmm::least_squares_cg(AM, UM, FM, iter);
		}
		// Check
		if (iter.converged())
			cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
		else if (iter.get_iteration() == descr_transp.MAXITER)
			cerr << "  ... reached the maximum number of iterations!" << endl;

	}*/
	
	if(NANO==1 && SATURATION==0)
	{
				
			//Compute Psi
			if(t==0){
			vector_type Cvold_proj(mf_coefv.nb_dof());gmm::clear(Cvold_proj);
			gmm::copy(UM_transp_v, Cvold);
			vector_type f_Cvold_proj(mf_coefv.nb_dof()); gmm::clear(f_Cvold_proj);
			

		    //Projection of Cvold on PK(1,0) for computing psi	
			//Build the rhs terms for the projection
			getfem::generic_assembly assem;
			assem.push_mi(mimv);
			assem.push_mf(mf_Cv);
			assem.push_mf(mf_coefv);
			assem.push_data(Cvold);
			assem.push_vec(f_Cvold_proj); 		
			assem.set("m1=data$1(#1);"
				  "V$1(#2)+=comp(Base(#1).Base(#2))(j,:).m1(j);");
			assem.assembly();

			gmm::csc_matrix<scalar_type> Mcc_aux1;
			gmm::clean(Mcc, 1E-12);
			gmm::copy(Mcc, Mcc_aux1);
			scalar_type cond;
			//solve the projection

			gmm::SuperLU_solve(Mcc_aux1, Cvold_proj, f_Cvold_proj,cond);
		
		}
		else{
		

		gmm::clear(Cvnew);

		gmm::copy(UM_transp_v, Cvnew);

		vector_type f_Cvnew_proj(mf_coefv.nb_dof()); gmm::clear(f_Cvnew_proj);
		vector_type Cvnew_proj(mf_coefv.nb_dof()); gmm::clear(Cvnew_proj);
		vector_type f_Cvold_proj(mf_coefv.nb_dof()); gmm::clear(f_Cvold_proj);
		vector_type Cvold_proj(mf_coefv.nb_dof());gmm::clear(Cvold_proj);

	    //Projection of Cvnew on PK(1,0) for computing psi	
		//Build the rhs terms for the projection
		getfem::generic_assembly assem;
		assem.push_mi(mimv);
		assem.push_mf(mf_Cv);
		assem.push_mf(mf_coefv);
		assem.push_data(Cvnew);
		assem.push_data(Cvold);
		assem.push_vec(f_Cvnew_proj); 	
		assem.push_vec(f_Cvold_proj);	
		assem.set("m1=data$1(#1);"
			  "V$1(#2)+=comp(Base(#1).Base(#2))(j,:).m1(j);");
		assem.set("m2=data$2(#1);"
			  "V$2(#2)+=comp(Base(#1).Base(#2))(j,:).m2(j);");
		assem.assembly();

		gmm::csc_matrix<scalar_type> Mcc_aux1;
		gmm::clean(Mcc, 1E-12);
		gmm::copy(Mcc, Mcc_aux1);
		scalar_type cond;
		//solve the projection
		gmm::SuperLU_solve(Mcc_aux1, Cvnew_proj, f_Cvnew_proj,cond);
		gmm::SuperLU_solve(Mcc_aux1, Cvold_proj, f_Cvold_proj,cond);

    		vector_type temp(mf_coefv.nb_dof()); gmm::clear(temp);
    		scalar_type cost;
    		cost = param_transp.dt()*param_transp.dadim()/(param_transp.Uadim()*2);   //dimensional time
    		gmm::copy(Psi, Psiold);

    		gmm::copy(gmm::scaled(Cvnew_proj, cost), temp);
    		gmm::add(gmm::scaled(Cvold_proj, cost), temp);

    		for(size_type i=0; i<mf_coefv.nb_dof(); i++) { //Psi is dimensional [1/m^2]
      		Psi[i] = PiGrecoproj[i]*param_transp.Uadim()*temp[i];
	
      		}

 		gmm::add(Psiold, Psi);
		gmm::copy(Cvnew_proj, Cvold_proj);
		gmm::copy(Cvnew,Cvold);
		}
	}
	
/*	if (NANO==1 && SATURATION==1){
		
		//Compute Psi
		gmm::copy(UM_transp_v, Cvold);
		for(size_type i=0; i<dof_transp.Cv(); i++) { 
      			Psi[i] = PiGrecoproj_cv[i]*std::max(0.0,(psi_max_vec[i]-Psiold[i]))*Cvold[i];
      		}
		gmm::scale(Psi, param_transp.dadim()*param_transp.dt()/param_transp.psi_max());
		gmm::add(Psiold, Psi);
		for(size_type i=0; i<dof_transp.Cv(); i++) { 
      			Pigrecostar[i] = PiGrecoproj_scaled_cv[i]*std::max(0.0,(psi_max_vec[i]-Psi[i]));
      		}
		gmm::scale(Pigrecostar, 1.0/param_transp.psi_max());
		
		gmm::clear(Psiold);
		gmm::add(Psi, Psiold);
	}
*/	

	if (NANO==1 && SATURATION==1){
		gmm::copy(UM_transp_v, Cvold);
	}

	#ifdef M3D1D_VERBOSE_
  	std::cout<<"Solving tissue problem"<<std::endl;
	#endif

	vector_type J(dof_transp.Cv()); gmm::clear(J);

//=======================================================================
	//Projecting Psi on Cv perchè J è di tipo mf_Cv
		std::cout<< "===================Start Proiezione Psi "<<std::endl;
		
		vector_type f_Psi_proj(mf_Cv.nb_dof()); gmm::clear(f_Psi_proj);
		vector_type Psi_proj(mf_Cv.nb_dof()); gmm::clear(Psi_proj);

		getfem::generic_assembly assem;
		assem.push_mi(mimv);
		assem.push_mf(mf_Cv);
		assem.push_mf(mf_coefv);
		assem.push_data(Psi);
		assem.push_vec(f_Psi_proj); 		
		assem.set("m1=data$1(#2);"
			  "V$1(#1)+=comp(Base(#2).Base(#1))(j,:).m1(j);");
		assem.assembly();

		sparse_matrix_type Mcc_cv(dof_transp.Cv(), dof_transp.Cv()); gmm::clear(Mcc_cv);
		getfem::asm_mass_matrix(Mcc_cv, mimv, mf_Cv);
		scalar_type cond;
		//solve the projection
		gmm::SuperLU_solve(Mcc_cv, Psi_proj, f_Psi_proj,cond);
				std::cout<< "===================End Proiezione Psi "<<std::endl;
//========================================================================

	//gmm::copy(Psi_proj, J);
	//Scalo Psi per la concentrazione di riferimento			
	gmm::copy(gmm::scaled(Psi_proj,2.39E14), J);

	scalar_type J_coeff;
	if(t==0)
		J_coeff=0.0;
	else{
		//J_coeff= param_transp.ups()*param_transp.b()*pow(t,param_transp.b()-1.0)*param_transp.cNP()*param_transp.VNP()/pow((pow(t,param_transp.b())+param_transp.ups()),2.0);
		//Fede: I think J shoud be adimensionalized w.r.t. cNP because then we export Ct adimensionalized w.r.t. cNP		
		J_coeff= (1/param_transp.dadim())*param_transp.ups()*param_transp.b()*pow(t,param_transp.b()-1.0)*param_transp.VNP()/pow((pow(t,param_transp.b())+param_transp.ups()),2.0);	
	}
	
	gmm::scale(J, J_coeff);	
	vector_type RHS_nano_(dof_transp.Ct()); gmm::clear(RHS_nano_);
	gmm::mult(Btv1,J,RHS_nano_);


	update_tissue(RHS_nano_);
	
	gmm::csc_matrix<scalar_type> A_transp_t;
	gmm::clean(AM_transp_t, 1E-12);
	gmm::copy(AM_temp_t, A_transp_t);
	
	vector_type F_transp_t(gmm::vect_size(FM_transp_t));
	gmm::clean(FM_transp_t, 1E-12);
	gmm::copy(FM_temp_t, F_transp_t);
	

	if ( descr_transp.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method to tissue problem... " << endl;
		#endif
		scalar_type cond_t;
		gmm::SuperLU_solve(A_transp_t, UM_transp_t, F_transp_t, cond_t);
		cout << "  Condition number (transport tissue problem): " << cond_t << endl;
	}


	//export solution
	std::cout<<" Solved! going to export..."<<std::endl;
	string time_suff = "";
	std::ostringstream convert;
	convert << time_count;
	time_suff = convert.str();

	//Compute the total c_t in the tissue
	outMeanCt<<mean_ct()<<std::endl;
	
	//Export the first iteration
	if(t==0){
	export_vtk(time_suff); 
	}
	
	//Export each 60 time steps
	//if(iteraz%10==0){ 
	export_vtk(time_suff);
	std::cout<<"exported! now new iteration..."<<std::endl;
	//}

	} //end of cycle over time 

	outMeanCt.close();
	
	export_vtk_nano();

	cout << endl<<"... time to solve : " << gmm::uclock_sec() - time << " seconds\n";

	return true;
 }; // end of solve
	
	
 void transport3d1d::export_vtk (const string & time_suff,const string & suff)
 {
  if (PARAM.int_value("VTK_EXPORT"))
  {
	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_transp.Ct()); 
	// Array of unknown dof of the network velocity
	vector_type Cv(dof_transp.Cv()); 

	//Copy solution
	gmm::copy(UM_transp_t, Ct);
	gmm::copy(UM_transp_v, Cv);


	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ct: adimensional with respect to cNP ..." << endl;
	#endif
	vtk_export exp_Ct(descr_transp.OUTPUT+"Ct"+suff+"_t"+time_suff+".vtk");
	exp_Ct.exporting(mf_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct, Ct, "Ct_adim");



	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv..." << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+"Cv"+suff+"_t"+time_suff+".vtk");
	exp_Cv.exporting(mf_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_Cv, Cv, "Cv");
 	

	
	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif
  }
 }; // end of export
 
   void transport3d1d::export_vtk_nano (const string & suff)
 {
  if (PARAM.int_value("VTK_EXPORT"))
  {
	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the fluidynamcal quantities to " << descr_transp.OUTPUT << " ..." << endl;
	#endif

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Reynolds_v ..." << endl;
	#endif

	vtk_export exp_Rev(descr_transp.OUTPUT+"Reproj.vtk");
	exp_Rev.exporting(mf_coefv);
	exp_Rev.write_mesh();
	exp_Rev.write_point_data(mf_coefv, Reproj, "Re"); 

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Probability of adhesion_v..." << endl;
	#endif

	vtk_export exp_padhv(descr_transp.OUTPUT+"Padhproj.vtk");
	exp_padhv.exporting(mf_coefv);
	exp_padhv.write_mesh();
	exp_padhv.write_point_data(mf_coefv, P_adhproj, "Padhv"); 
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Wall shear stress_v ..." << endl;
	#endif
		
	vtk_export exp_WSRv(descr_transp.OUTPUT+"WSRproj.vtk");
	exp_WSRv.exporting(mf_coefv);
	exp_WSRv.write_mesh();
	exp_WSRv.write_point_data(mf_coefv, WSRproj, "WSRv"); 
		
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting vascular adhesion parameter Pi ..." << endl;
	#endif
		
	vtk_export exp_pigrecov(descr_transp.OUTPUT+"PiGrecoproj.vtk");
	exp_pigrecov.exporting(mf_coefv);
	exp_pigrecov.write_mesh();
	exp_pigrecov.write_point_data(mf_coefv, PiGrecoproj, "PiGreco_adim"); 

		
	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif
  }
 };
  
	
 } // end of namespace
