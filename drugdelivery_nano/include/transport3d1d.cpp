/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   transport3d1d.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 - May 2018
  @brief  Definition of the main class for the 3D/1D coupled transport problem.
 */
 
 #include <transport3d1d.hpp>
 #include <AMG_Interface.hpp>
 #include <cmath>
 #include "gmm/gmm_inoutput.h"
 #include "getfem/getfem_import.h"

#define M3D1D_VERBOSE_

 namespace getfem {


	// Initialize the transport problem
 	void transport3d1d::init_transp (int argc, char *argv[]) 
 	{
 	#ifdef M3D1D_VERBOSE_
 	std::cout << "initialize transport problem..."<<std::endl<<std::endl;
 	#endif

	PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	import_data_transp();
	//2. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh_transp();
	//3. Set finite elements and integration methods
	set_im_and_fem_transp();
 	//4. Build problem parameters
 	build_param_transp();
	//5. Build the list of tissue boundary data
	build_tissue_boundary_transp();
	//6. Build the list of tissue boundary (and junction) data
 	build_vessel_boundary_transp();

 	}; // end of init


 	// Aux methods for init
	
	// Import algorithm specifications
	void transport3d1d::import_data_transp(void)
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Importing descriptors for tissue and vessel problems ..." << endl;
	#endif

	descr.import(PARAM);
	descr_transp.import(PARAM);
	#ifdef M3D1D_VERBOSE_
	cout << descr_transp;
	#endif
	 
 
	}; //end of import_data_transp
	
	 
	
	// Import mesh for tissue (3D) and vessel (1D)  
	void transport3d1d::build_mesh_transp(void){

	//In order to have the boundary conditions for the nodes
	//we need to build again the 1D mesh from another pts file
	//! \todo write import_msh_file_transp which take in account both .pts files: modify transport3d1d::init_fluid such that it calls import_msh_file_transp.
	//! \todo branch index should start from 1, not from 0: the junctions will use the notation that +branch is inflow and -branch is outflow. Branch 0 force the user to not give in the .pts file an outflow for first branch (this is usually ok, but loses generality)
	mesht.clear();
		bool test = 0;
	test = PARAM.int_value("TEST_GEOMETRY");
	if(test==0){
		#ifdef M3D1D_VERBOSE_
		cout << "Importing the 3D mesh for the tissue ...  "   << endl;
		#endif
	string st("gmsh:"+descr.MESH_FILET);
	getfem::import_mesh(st,mesht);
		 //import_msh_file(descr.MESH_FILET, mesht);
	}else{
		#ifdef M3D1D_VERBOSE_
		cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
		#endif
		string st("GT='" + PARAM.string_value("GT_T") + "'; " +
					   "NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
					   "ORG=" + PARAM.string_value("ORG_T") + "; " +  
					   "SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
					   "NOISED=" + PARAM.string_value("NOISED_T")); 
		#ifdef M3D1D_VERBOSE_		
		cout << "mesht description: " << st << endl;
		#endif
		regular_mesh(problem3d1d::mesht, st);
	}
	
	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel (transport problem)... "   << endl;
	#endif
	std::ifstream ifs(descr_transp.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr_transp.MESH_FILEV);
	meshv.clear();
	bool Import=PARAM.int_value("IMPORT_CURVE");
	bool Curve=PARAM.int_value("CURVE_PROBLEM");

	if(Curve && !Import){
		import_pts_file(ifs, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV, param);
	}
	else if(Import && !Curve){
		GMM_ASSERT1(0,"If you want to import the curvature, you need to enable CURVE_PROBLEM=1");
	}
	else if(Import && Curve){
		std::ifstream ifc(PARAM.string_value("CURVE_FILE","curvature file location"));
		GMM_ASSERT1(ifc.good(), "impossible to read from file " << PARAM.string_value("CURVE_FILE","curvature file location"));
		
		import_pts_file(ifs,ifc, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV, param);

		ifc.close();
	} else{
		import_pts_file(ifs, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV);
	}


	nb_branches = nb_vertices.size();
	ifs.close();

	
	}; // end of build_mesh_geometry

	// Set finite elements methods and integration methods 
	void transport3d1d::set_im_and_fem_transp(void)
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
	

	// I had to delete the meshes and build them again, because of issues on the boundary conditions:
	// Now i have to build again also the mesh_fem and mesh_im objects.
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

	}; // end of set_im_and_fem
	
	
	// Build problem parameters
	void transport3d1d::build_param_transp(void)
	{
	
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif
	param.build(PARAM, mf_coeft, mf_coefv,mf_coefvi); 
	param_transp.build(PARAM, mf_coeft, mf_coefv);
	#ifdef M3D1D_VERBOSE_
	cout << param_transp ;
	#endif

	}; // end of build_param_transp
  
  
  	//Build boundary regions on tissue
	void
	transport3d1d::build_tissue_boundary_transp (void) 
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Building tissue boundary ..." << endl;
	#endif

	size_type face=descr_transp.FACE;

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
		BCt_transp.emplace_back(labels[f], std::stof(values[f]), 0, face+f);
		#ifdef M3D1D_VERBOSE_
		cout << "  face " << f << " : " << BCt_transp.back() << endl;
		#endif
	} 
	

	// Build mesht regions
	size_type xx=0, yy=1, zz=2;

	mesh_region border_faces;
	outer_faces_of_mesh(mesht, border_faces);

	for (mr_visitor i(border_faces); !i.finished(); ++i) {

		assert(i.is_face());

		// Unit outward normal : used to identify faces
		//! \todo Use getfem 5.0's function select_faces_of_normal?
		base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		     if (gmm::abs(un[xx] + 1.0) < 1.0E-7) 	// back
			mesht.region(face+0).add(i.cv(), i.f());
		else if (gmm::abs(un[xx] - 1.0) < 1.0E-7) 	// front
			mesht.region(face+1).add(i.cv(), i.f());
		else if (gmm::abs(un[yy] + 1.0) < 1.0E-7) 	// left
			mesht.region(face+2).add(i.cv(), i.f());
		else if (gmm::abs(un[yy] - 1.0) < 1.0E-7) 	// right
			mesht.region(face+3).add(i.cv(), i.f());
		else if (gmm::abs(un[zz] + 1.0) < 1.0E-7) 	// bottom
			mesht.region(face+4).add(i.cv(), i.f());
		else if (gmm::abs(un[zz] - 1.0) < 1.0E-7) 	// top
			mesht.region(face+5).add(i.cv(), i.f());
				
		}



	GMM_ASSERT1(!(PARAM.int_value("TEST_GEOMETRY") && descr_transp.CONFORMING), "Test 3D mesh cannot be conforming! Check the .param file!!!");
	
	if(descr_transp.CONFORMING){
	/*! \todo Regions SIGMA and GAMMA are described in .msh file from GMESH (using PhysicalVolume). 
                  Build directly those regions ("if(dist(element, meshv)>0){OMEGA.add(element)}") 
		  or build the indicator functions and rescale all the integrals by those functions.
	*/
	#ifdef M3D1D_VERBOSE_
	cout << "The mesh is conforming with respect to the vessel!" << endl;
	cout << "Checking the regions..." << endl;
	#endif

	size_type sigma=descr_transp.SIGMA;
	size_type omega=descr_transp.OMEGA;
	size_type gamma=descr_transp.GAMMA;	

	//check if there is a overlap in regions!
	GMM_ASSERT1(sigma!=gamma, "SIGMA=GAMMA: check the .param file");
	GMM_ASSERT1(omega!=sigma, "OMEGA=SIGMA: check the .param file");
	GMM_ASSERT1(gamma!=omega, "GAMMA=OMEGA: check the .param file");
	
	for(int i=0; i<5; i++){
	GMM_ASSERT1((face+i)!=gamma, "FACE=GAMMA: check the .param file");
	GMM_ASSERT1((face+i)!=sigma, "FACE=SIGMA: check the .param file");
	GMM_ASSERT1((face+i)!=omega, "FACE=OMEGA: check the .param file");
}

	//Check if sigma and omega are defined in the msh file
	GMM_ASSERT1(mesht.has_region(sigma), "File .msh does not contain region SIGMA! Check .msh and .param!!");
	GMM_ASSERT1(mesht.has_region(omega), "File .msh does not contain region OMEGA! Check .msh and .param!!");

	// GAMMA is the boundary of the vessel SIGMA:
	
	mesh_region gamma_region;
	outer_faces_of_mesh(mesht, mesht.region(omega), gamma_region);

	for (mr_visitor i(gamma_region); !i.finished(); ++i) {

		assert(i.is_face());

		     if (!(mesht.region(face+0).is_in(i.cv(),i.f()))&&!(mesht.region(face+1).is_in(i.cv(),i.f()))
				&&!(mesht.region(face+2).is_in(i.cv(),i.f()))&&!(mesht.region(face+3).is_in(i.cv(),i.f()))
				&&!(mesht.region(face+4).is_in(i.cv(),i.f()))&&!(mesht.region(face+5).is_in(i.cv(),i.f()))) 	// back
			mesht.region(gamma).add(i.cv(), i.f());
			mesht.region(gamma+1).add(i.cv(), i.f());
				
		} 

	mesh_region gamma_region2;
	outer_faces_of_mesh(mesht, mesht.region(sigma), gamma_region2);

	for (mr_visitor i(gamma_region2); !i.finished(); ++i) {

		assert(i.is_face());

		     if (!(mesht.region(face+0).is_in(i.cv(),i.f()))&&!(mesht.region(face+1).is_in(i.cv(),i.f()))
				&&!(mesht.region(face+2).is_in(i.cv(),i.f()))&&!(mesht.region(face+3).is_in(i.cv(),i.f()))
				&&!(mesht.region(face+4).is_in(i.cv(),i.f()))&&!(mesht.region(face+5).is_in(i.cv(),i.f()))) 	// back
			mesht.region(gamma).add(i.cv(), i.f());	
			mesht.region(gamma+2).add(i.cv(), i.f());
			//mesht.region(gamma+1).add(i.cv(), i.f());					
		}

	//outer_faces_of_mesh(mesht, mesht.region(sigma), mesht.region(gamma));
	GMM_ASSERT1(mesht.has_region(gamma), "Something wrong happened: I couldn't build the region GAMMA");

	outer_faces_of_mesh(mesht, mesht.region(omega), mesht.region(gamma+1));
	GMM_ASSERT1(mesht.has_region(gamma), "Something wrong happened: I couldn't build the region GAMMA+1");


	#ifdef M3D1D_VERBOSE_
	cout << "...Check complete! All 3D regions are correctly defined!" << endl;
	#endif		

}



	}; //end of build_tissue_boundary_transp

	//Build boundary regions on network

	void 
	transport3d1d::build_vessel_boundary_transp(void)
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
	#endif
	try {

	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector extrema;   // global idx of extreme vertices in meshv

	Jv_transp.clear();
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
		// Check if junction has been already stored, 
			// if not add to the junction list (J) and build a new region

	//! \todo add multiple times the junction node to the junction region. Generic_assembly looks (apparently) only at indexes of the nodes, not at his coordinates; in this way, when I build the region with the junction node from a certain branch, generic_assembly will not recognize the same node from another branch (probably i look at the basis buildt only on the first branch). In order to use the generic_assembly for junction nodes I should add all the basis to the region (e.g. the same node from all the branches)

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
				Jv_transp.emplace_back("JUN", 0, i0, fer);
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
				found = (i0 == Jv_transp[jj].idx);
				if (!found) jj++;
			}
			//cout << "Branch -" << branch << " added to junction " << jj << endl;
			Jv_transp[jj].value += param.R(mimv, branch);
			Jv_transp[jj].branches.emplace_back(-branch);
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

			else { // interior -> Mixed point 
				// "MIX" label via post-processing
			// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				BCv_transp.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_transp.back().branches.emplace_back(branch); 
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
					Jv_transp.emplace_back("JUN", 0, i1, fer);
					fer++;
				// Search for index of second containing branch (\mathcal{P}^{out}_j)
				size_type secondbranch = 0; 
				size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
				size_type firstcv = (( cv1 != cv) ? cv2 : cv1);
				contained = false;
				while (!contained && secondbranch<nb_branches ) {
					if (secondbranch!=firstbranch)
					contained = meshv.region(secondbranch).is_in(secondcv);
					if (!contained) secondbranch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				// Add the two branches
				scalar_type in;
				in=0;
				if (meshv.ind_points_of_convex(firstcv)[0]==i1) in=-1;
				else if (meshv.ind_points_of_convex(firstcv)[1]==i1) in=+1;
				GMM_ASSERT1(in!=0, "There's something wrong in firstbranch convex index");
				Jv_transp.back().branches.emplace_back(in*firstbranch);

				in=0;
				if (meshv.ind_points_of_convex(secondcv)[0]==i1) in=-1;
				else if (meshv.ind_points_of_convex(secondcv)[1]==i1) in=+1;
				GMM_ASSERT1(in!=0, "There's something wrong in secondbranch convex index");
				Jv_transp.back().branches.emplace_back(in*secondbranch);
				Jv_transp.back().value += param.R(mimv, firstbranch);
				Jv_transp.back().value += param.R(mimv, secondbranch);
				}
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
				Jv_transp.emplace_back("JUN", 0, i1, fer);
				// Add the outflow branch
				Jv_transp.back().branches.emplace_back(+branch);
				Jv_transp.back().value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << i1 << endl;
				fer++;
			}
			else {
				// Add the outflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i1 == Jv_transp[jj].idx);
					if (!found) jj++;
				}
				Jv_transp[jj].branches.emplace_back(+branch);
				Jv_transp[jj].value += param.R(mimv, branch);
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
	for (size_type i=0; i<Jv_transp.size(); ++i)
		cout << "    -  label=" << Jv_transp[i].label 
			 << ", value=" << Jv_transp[i].value << ", ind=" << Jv_transp[i].idx 
			 << ", rg=" << Jv_transp[i].rg << ", branches=" << Jv_transp[i].branches << endl; 
	cout << "---------------------------------------- "   << endl;
	#endif

	} 
	GMM_STANDARD_CATCH_ERROR; // catches standard errors

	} /* end of build_vessel_boundary_transp */


  
	void transport3d1d::assembly_transp (void)
	 {

 	 //Build the monolithic matrix AM
	 assembly_mat_tissue_transp();
	 assembly_mat_network_transp();
	 
	 }; // end of assembly

	void transport3d1d::assembly_mat_tissue_transp(void)
		{
		#ifdef M3D1D_VERBOSE_
		cout << "Allocating AM_transp, UM_transp, FM_transp ..." << endl;
		#endif
		gmm::resize(AM_transp_t, dof_transp.Ct(), dof_transp.Ct());	gmm::clear(AM_transp_t);
		gmm::resize(UM_transp_t, dof_transp.Ct()); 	gmm::clear(UM_transp_t);
		gmm::resize(FM_transp_t, dof_transp.Ct()); 	gmm::clear(FM_transp_t);
		
		
		#ifdef M3D1D_VERBOSE_
		cout << "Assembling the matrix AM_transp_t for tissue ..." << endl;
		#endif
		// Mass(time derivative)  matrix for the interstitial problem
		sparse_matrix_type Mt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Mt);
		// Diffusion matrix for the interstitial problem
		sparse_matrix_type Dt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Dt);
		//Advection matrix for interstitial problem
		sparse_matrix_type Bt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Bt);
		// Reaction matrix for the interstitial problem
		sparse_matrix_type Rt(dof_transp.Ct(), dof_transp.Ct()); gmm::clear(Rt);

		#ifdef M3D1D_VERBOSE_
		cout<< "   Computing Peclet number..."<<endl;	
		#endif	
		//vectors containing the exact solution
		vector_type Ut(dof.Ut()); 	
		gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())), Ut);
		//Interpolate Ut on polinomials of degree 0 
		pfem pf_U = fem_descriptor("FEM_PK(3,0)");
		mesh_fem mf_U(mesht);
		mf_U.set_qdim(3);
		mf_U.set_finite_element(mesht.convex_index(), pf_U);
		vector_type Ut_(mf_U.nb_dof());
		getfem::interpolation(mf_Ut, mf_U, Ut, Ut_);
		//compute peclet
		scalar_type peclet_t= peclet(mesht, Ut_, param_transp.At(1), 3);

		cout<< "Peclet in tissue:    "<< peclet_t<<endl;

		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling Mt, Dt, Rt ..." << endl;
		#endif
		
		//Coefficient for mass term:
		vector_type mass_coeff(dof.Pt()); gmm::clear(mass_coeff);
		vector_type Pl(dof.Pt(),PARAM.real_value("PL"));
		gmm::scale(Pl, -1.0); 
		gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut(), dof.Pt())) ,  mass_coeff);
		gmm::add(Pl ,  mass_coeff);
		gmm::scale (mass_coeff,param_transp.Q_pl(1));

		//Build Mt, Dt and Rt
		asm_tissue_transp(Mt, Dt, Rt, mimt, mf_Ct, mf_coeft, param_transp.At(), mass_coeff);

		// Copy Mt: Time Derivative in tissue
		if(descr_transp.STATIONARY ==0)
		{ 
			gmm::scale(Mt, (1.0/param_transp.dt()));
			gmm::add(Mt, AM_transp_t); 
		}
			
		// Check peclet number for instability
		if((descr_transp.ADVECTION==1) && (peclet_t>1))
		{ 	
			cout<<"WARNING!! Peclet > 1 in tissue: applying artificial diffusion"<<std::endl;	
		  	gmm::scale(Dt, (1+peclet_t));
		}
		
		// Copy Dt: diffusion in tissue		  
		gmm::add(Dt, AM_transp_t);  	
		// Copy Rt: reaction in tissue
	 	gmm::add(Rt, AM_transp_t); 

		
		if(descr_transp.ADVECTION ==0)	
			{cout<<"No advection: only diffusion and reaction terms"<<endl;}
		else{			
			#ifdef M3D1D_VERBOSE_
			cout << "  Assembling Bt ..." << endl;
			#endif			
			// Build Bt
			asm_advection_tissue(Bt, mimt, mf_Ct, mf_Ut, gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())));
			// Copy Bt: advection in tissue
			gmm::add(Bt, AM_transp_t); 	
		}

		
		#ifdef M3D1D_VERBOSE_
		cout << "  Setting initial condition for tissue and network concentration ..." << endl;
		#endif

		vector_type C0t_vect(dof_transp.Ct(), param_transp.C0t());
		gmm::copy(C0t_vect, UM_transp_t);

		gmm::clear(C0t_vect);	

		// De-allocate memory
		gmm::clear(Mt);    		
		gmm::clear(Dt);  
		gmm::clear(Bt);    
		
	} /* end of assembly_mat_transp */

	void transport3d1d::assembly_mat_network_transp()
	{
		#ifdef M3D1D_VERBOSE_
			cout << "Allocating AM_transp, UM_transp, FM_transp ..." << endl;
		#endif

		gmm::resize(AM_transp_v, dof_transp.Cv(), dof_transp.Cv());	gmm::clear(AM_transp_v);
		gmm::resize(UM_transp_v, dof_transp.Cv()); 	gmm::clear(UM_transp_v);
		gmm::resize(FM_transp_v, dof_transp.Cv()); 	gmm::clear(FM_transp_v);
		
		#ifdef M3D1D_VERBOSE_
		cout << "Assembling the monolithic matrix AM_transp_v ..." << endl;
		#endif
		
		// Mass (time derivative)  matrix for the network problem
		sparse_matrix_type Mv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mv);		
		// Diffusion matrix for the network problem
		sparse_matrix_type Dv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Dv);
		//Advection matrix for network problem
		sparse_matrix_type Bv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bv);
		// Adhesion  matrix for the network problem					
		sparse_matrix_type Adhv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Adhv);	

		#ifdef M3D1D_VERBOSE_
			cout<< "   Computing Peclet number..."<<endl;	
		#endif	
		//vectors containing the exact solution
		vector_type Uv(dof.Uv());
		gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())), Uv);
		//compute peclet
		scalar_type peclet_v= peclet(meshv, Uv, param_transp.Av(1), 1);
		
		#ifdef M3D1D_VERBOSE_
		cout<< "Peclet in vessels:    "<< peclet_v<<endl;
		#endif	
				
		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling Mv and Dv ..." << endl;
		#endif	

		// Build Mv and Dv
		asm_network_transp(Mv, Dv, mimv,mf_Cv, mf_coefv, param_transp.Av(), param.R());
			
		// Copy Mv: Time Derivative in network
		if(descr_transp.STATIONARY ==0)
		{
			gmm::scale(Mv, (1.0/param_transp.dt()));
			gmm::add(Mv, AM_transp_v); 
		}

		// Check peclet number for instability
		if((descr_transp.ADVECTION==1) && (peclet_v>1))
			{ cout<<"WARNING!! Peclet > 1 in network: applying artificial diffusion"<<endl;
	   	 	  gmm::scale(Dv, (1+peclet_v)); }
			
		// Copy Dv: diffusion in network		 	
		gmm::add(Dv, AM_transp_v);
		
		if(descr_transp.ADVECTION ==0)	{cout<<"No advection: only diffusion and reaction terms"<<endl;}
		else
		{			
			#ifdef M3D1D_VERBOSE_
				cout << "  Assembling Bv ..." << endl;
			#endif	

			size_type shift = 0;
			for(size_type i=0; i<nb_branches; ++i)
			{
				if(i>0) shift += mf_Uvi[i-1].nb_dof();
				vector_type Uvi( mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
				gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);
				asm_advection_network(Bv, mimv, mf_Cv, mf_coefvi[i], mf_Uvi[i], mf_coefv, Uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i), param.R(), meshv.region(i));
			}
			gmm::scale(Bv, pi);
			// Copy Bv: advection in network
			gmm::add(Bv, AM_transp_v);
		
		}

		// Build the adhesion term 

		///////////////////////////////////////
		//Projection
		//////////////////////////////////////
		//Create the vectors for the projection on the coefv_mesh (global on the network)
		gmm::resize(Reproj, mf_coefv.nb_dof()); gmm::clear(Reproj);
		gmm::resize(WSRproj, mf_coefv.nb_dof()); gmm::clear(WSRproj);
		gmm::resize(P_adhproj, mf_coefv.nb_dof()); gmm::clear(P_adhproj);
		gmm::resize(PiGrecoproj, mf_coefv.nb_dof()); gmm::clear(PiGrecoproj);
		gmm::resize(PiGrecoproj_scaled, mf_coefv.nb_dof()); gmm::clear(PiGrecoproj_scaled);	
		gmm::resize(Uv_abs, mf_coefv.nb_dof()); gmm::clear(Uv_abs);	


		//mass matrix for the coefficient
		gmm::resize (Mcc,mf_coefv.nb_dof(), mf_coefv.nb_dof()); gmm::clear(Mcc);
		getfem::asm_mass_matrix(Mcc, mimv, mf_coefv);

		//rhs term
		vector_type f_Re(mf_coefv.nb_dof()); gmm::clear(f_Re);
		vector_type f_WSR(mf_coefv.nb_dof()); gmm::clear(f_WSR);
		vector_type f_Padh(mf_coefv.nb_dof()); gmm::clear(f_Padh);
		vector_type f_Pigreco(mf_coefv.nb_dof()); gmm::clear(f_Pigreco);
		vector_type f_Pigreco_scaled(mf_coefv.nb_dof()); gmm::clear(f_Pigreco_scaled);
		vector_type f_velocity(mf_coefv.nb_dof()); gmm::clear(f_velocity);
		
		bool NANO = PARAM.int_value("NANO", "flag for adhesive term in vessel");
		bool SATURATION  = PARAM.int_value("SATURATION","Flag to use the saturation per the adhesive term");

		size_type shift1 = 0; //dof of mf_Uvi[i]
		size_type shift2 = 0; //dof of mf_coefvi[i]

		for (size_type i=0; i<nb_branches; ++i){
			// Wall shear rate vector on a branch
			vector_type WSRi (mf_coefvi[i].nb_dof()); gmm::clear(WSRi);
			// Probability of adhesion vector on a branch
			vector_type P_adhi (mf_coefvi[i].nb_dof()); gmm::clear(P_adhi); 
			// Reynolds vector on a branch
			vector_type Rei (mf_coefvi[i].nb_dof()); gmm::clear(Rei);
		
			if(i>0) shift1 += mf_Uvi[i-1].nb_dof();
			if(i>0) shift2 += mf_coefvi[i-1].nb_dof();
			// Velocity of the i-th branch
	 		vector_type Uvi1(mf_Uvi[i].nb_dof()); gmm::clear(Uvi1);
			vector_type Uvi1_abs(mf_Uvi[i].nb_dof()); gmm::clear(Uvi1_abs);

			gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift1, mf_Uvi[i].nb_dof())) ,  Uvi1);


			for (size_type k=0; k<mf_Uvi[i].nb_dof(); ++k){
				Uvi1_abs[k]=fabs(Uvi1[k]);
			}

			//Reynolds and WSR should be defined on mf_coefv as the viscosity from problemHT
			//for this reason we project the velocity on mf_coefvi
			vector_type Uvi1abs_proj (mf_coefvi[i].nb_dof()); gmm::clear(Uvi1abs_proj);
			vector_type f_Uvi1 (mf_coefvi[i].nb_dof()); gmm::clear(f_Uvi1);
			getfem::generic_assembly assemm;
			assemm.push_mi(mimv);
			assemm.push_mf(mf_Uvi[i]);
			assemm.push_mf(mf_coefvi[i]);
			assemm.push_data(Uvi1_abs);
			assemm.push_vec(f_Uvi1);
			assemm.set("m1=data$1(#1);"
				  "V$1(#2)+=comp(Base(#1).Base(#2))(j,:).m1(j);"
				); 
			assemm.assembly();	
				
			sparse_matrix_type Proj_temp (mf_coefvi[i].nb_dof(),mf_coefvi[i].nb_dof()); gmm::clear(Proj_temp);
			asm_mass_matrix(Proj_temp, mimv, mf_coefvi[i]);
			scalar_type condd;
			gmm::SuperLU_solve(Proj_temp, Uvi1abs_proj, f_Uvi1,condd);

			
			// Radius of the i-th branch
			scalar_type Ri = param.R(mimv, i); //raggio della regione (ramo) i
			
			// Wall shear rate		
			wall_shear_stress(WSRi, Ri, Uvi1abs_proj); //Adimensional for the export 
			vector_type WSRi_dim(mf_coefvi[i].nb_dof()); gmm::clear(WSRi_dim);  //to be put in P_a
			gmm::add(WSRi,WSRi_dim);	
			gmm::scale(WSRi_dim, (param_transp.Uadim()/ param_transp.dadim()));// Rescale WSRi_adim * U / d -> WSR is dimensional
			// Reynolds
			//get the viscosity on the branch i from problemHT
			vector_type mu_HTi(mf_coefvi[i].nb_dof()); gmm::clear(mu_HTi); // dimensional
			if(descrHT.HEMATOCRIT_TRANS){
			size_type pos=0;
				for (getfem::mr_visitor mrv(mf_coefv.linked_mesh().region(i)); !mrv.finished(); ++mrv)
					for (auto muu : mf_coefv.ind_basic_dof_of_element(mrv.cv()))
						{	
							mu_HTi[pos] = MU[muu];
							pos++;}
						}
			else{
				for(size_type k = 0; k< mf_coefvi[i].nb_dof(); k++)
					mu_HTi[k] = param.mu_v();
			}
			reynolds (Rei, Uvi1abs_proj, Ri, mu_HTi, param_transp.rho()); // Dimensional radius and velocity 
			gmm::scale(Rei,(param_transp.Uadim()*param_transp.dadim())); //Rescale Re_adimensi * U * d -> Rei is the real Reynolds
			//Probability of adhesion
			bool PA_INPUT  = PARAM.int_value("PA_INPUT");
			//If probability of adhesion from the exponential function
			if(PA_INPUT==0){
			scalar_type alpha2=param_transp.m_r()*(1.0-pow(1.0-(2.0*param_transp.h_0()/param_transp.dp()),2.0));
			probability_adhesion0(P_adhi, WSRi_dim, param_transp.m_l(), param_transp.Ka(),alpha2,param_transp.r_0(),param_transp.beta_nano(),mu_HTi,param_transp.m_r()); //Dimensional WRS
				}

			// If probability of adhesion from Lattice Boltzmann method
			// Case strong bond -> 
			//! \todo Modify this and pass these values as input
			if(PA_INPUT==1){
			scalar_type Pa_min;
			scalar_type Pa_med;
			scalar_type Pa_max;
			//------------------------------------------------------------------------------
					//Values of Pa for rigid strong bond
				if(param_transp.rho_l()==0.3){Pa_min=0.3; Pa_med=0.0; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.5){Pa_min=1.0; Pa_med=0.33; Pa_max=0.28;}
				else if (param_transp.rho_l()==0.7){Pa_min=1.0; Pa_med=0.64; Pa_max=0.56;}
				else if (param_transp.rho_l()==0.9){Pa_min=1.0; Pa_med=0.8; Pa_max=0.7;}
				else {GMM_ASSERT1(PA_INPUT==1, "rho_l value not valid");}
			//------------------------------------------------------------------------------
				//Values of Pa for rigid mild bond
				/*if(param_transp.rho_l()==0.3){Pa_min=0.0; Pa_med=0.0; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.5){Pa_min=0.56; Pa_med=0.05; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.7){Pa_min=1.0; Pa_med=0.64; Pa_max=0.56;}
				else if (param_transp.rho_l()==0.9){Pa_min=1.0; Pa_med=0.8; Pa_max=0.7;}
				else {GMM_ASSERT1(PA_INPUT==1, "rho_l value not valid");}*/
			//------------------------------------------------------------------------------
				//Values of Pa for soft strong bond
				/*if(param_transp.rho_l()==0.3){Pa_min=0.5349077; Pa_med=0.3747868; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.5){Pa_min=0.7594362; Pa_med=0.4155242; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.7){Pa_min=0.7886917; Pa_med=0.4292456; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.9){Pa_min=1.0; Pa_med=0.5895871; Pa_max=0.0;}
				else {GMM_ASSERT1(PA_INPUT==1, "rho_l value not valid");}*/
			//------------------------------------------------------------------------------
				//Values of Pa for soft mild bond
				/*if(param_transp.rho_l()==0.3){Pa_min=0.6881118; Pa_med=0.0; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.5){Pa_min=0.7068774; Pa_med=0.0; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.7){Pa_min=0.7786835; Pa_med=0.0; Pa_max=0.0;}
				else if (param_transp.rho_l()==0.9){Pa_min=1.0; Pa_med=0.0; Pa_max=0.0;}
				else {GMM_ASSERT1(PA_INPUT==1, "rho_l value not valid");}*/
			probability_adhesion(P_adhi, Rei, Pa_min, Pa_med, Pa_max);	
				}
			
			vector_type adhesion_coeff(mf_coefvi[i].nb_dof()); gmm::clear(adhesion_coeff); 
			vector_type adhesion_coeff_scaled(mf_coefvi[i].nb_dof()); gmm::clear(adhesion_coeff_scaled); 
			vector_type temporary(mf_coefvi[i].nb_dof()); gmm::clear(temporary); 
			gmm::add(WSRi_dim, temporary); 
			for (size_type i =0; i< P_adhi.size(); i++)				// P_adh * WSR (component by component) -> P_adh; 
				temporary[i] = P_adhi[i] * temporary[i]; 
			gmm::scale(temporary, (param_transp.dp()/2.0)); 	// Scale dp/2
			gmm::scale(temporary, (1.0/param_transp.Uadim())); 	// Adimensionalization
			gmm::add(temporary,adhesion_coeff);

			gmm::add(temporary,adhesion_coeff_scaled);
			//gmm::scale(adhesion_coeff_scaled, 2.0/Ri); 
			gmm::scale(adhesion_coeff_scaled, 2.0*pi*Ri);
			
			if(NANO==1 && SATURATION==0){
			asm_network_nano_transp(Adhv, mimv, mf_Cv, mf_coefvi[i], adhesion_coeff_scaled, meshv.region(i));
			masslumping(Adhv);	
				}

			//Build the rhs terms for the projection
			getfem::generic_assembly assem;
			assem.push_mi(mimv);
			//assem.push_mf(mf_Uvi[i]);
			assem.push_mf(mf_coefvi[i]);
			assem.push_mf(mf_coefv);
		     	// assem.push_mf(mf_Cv);
			assem.push_data(Rei);
			assem.push_data(WSRi);
			assem.push_data(P_adhi);
			assem.push_data(adhesion_coeff);
			assem.push_data(adhesion_coeff_scaled);
			assem.push_data(Uvi1abs_proj); //for velocity export
			assem.push_vec(f_Re);
			assem.push_vec(f_WSR);
			assem.push_vec(f_Padh);
			assem.push_vec(f_Pigreco);
			assem.push_vec(f_Pigreco_scaled);
			assem.push_vec(f_velocity);
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
				  "m6=data$6(#1);"
				  "V$6(#2)+=comp(Base(#1).Base(#2))(j,:).m6(j);"
				); 
			assem.assembly(meshv.region(i));

			
		}

		// solve the projection problem
		#ifdef M3D1D_VERBOSE_
		cout << "  Solve the projection of the fluid dynamical quantities ... " << endl;
		#endif
		gmm::csc_matrix<scalar_type> Mcc_aux;
		gmm::clean(Mcc, 1E-12);
		gmm::copy(Mcc, Mcc_aux);
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
		//Vascular adhesion parameter PiGreco*2/R (adimensional) on coefv--> 2/R*PiGreco/Uadim
		gmm::SuperLU_solve(Mcc_aux, Uv_abs, f_velocity,cond);
		
		gmm::clear(f_Re);
		gmm::clear(f_WSR);
		gmm::clear(f_Padh);
		gmm::clear(f_Pigreco);
		gmm::clear(f_Pigreco_scaled);
		gmm::clear(f_velocity);

		//Adding adhesion term
		if(NANO==1 && SATURATION==0)
			gmm::add(Adhv, AM_transp_v);
		

		#ifdef M3D1D_VERBOSE_
			cout << "  Setting initial condition for tissue and network concentration ..." << endl;
		#endif

		vector_type C0v_vect(dof_transp.Cv(), param_transp.C0v());
		gmm::copy(C0v_vect, UM_transp_v);
		gmm::clear(C0v_vect);

		// De-allocate memory
	    gmm::clear(Mv); 
		gmm::clear(Dv);
		gmm::clear(Bv);
		
	} /* end of assembly_mat_transp */

	void transport3d1d::assembly_rhs_tissue_transp()
	{
 
		#ifdef M3D1D_VERBOSE_
			cout << "Assembling the monolithic rhs FM_transp_t ... " << endl;
		#endif

		//Right Hand Side for tissue	
		vector_type Ft(dof_transp.Ct());
		sparse_matrix_type Att(dof_transp.Ct(), dof_transp.Ct());

		gmm::add(AM_temp_t, Att);
		gmm::scale(AM_temp_t, 0.0);	
				
		gmm::add(FM_temp_t,Ft);	 
		gmm::scale(FM_temp_t,0.0);
		
		scalar_type beta_t  = PARAM.real_value("BETAtissue_transp", "Coefficient for mixed BC for transport problem in tissue");
		asm_tissue_bc_transp(Ft, Att, mimt, mf_Ct, mf_coeft, BCt_transp,beta_t);
		gmm::add(Att,AM_temp_t);
		gmm::add(Ft, FM_temp_t);
		// De-allocate memory
		gmm::clear(Att);
		gmm::clear(Ft);
	}

	void transport3d1d::assembly_rhs_network_transp(){
		#ifdef M3D1D_VERBOSE_
			cout << "  Building vessel boundary term ..." << endl;
		#endif
			
		//Right Hand Side for vessels
		sparse_matrix_type Avv(dof_transp.Cv(), dof_transp.Cv());
		vector_type Fv(dof_transp.Cv());
		gmm::add(AM_temp_v, Avv);
		gmm::scale(AM_temp_v, 0.0);				
		gmm::add(FM_temp_v,Fv);	
		gmm::scale(FM_temp_v,0.0);

		scalar_type beta_v  = PARAM.real_value("BETAvessel_transp", "Coefficient for mixed BC for transport problem in vessels");
		if(t<=param_transp.inj_time())
			asm_network_bc_transp(Fv, Avv, mimv, mf_Cv, mf_coefv, BCv_transp, beta_v, param.R() );
		else{
			for(size_type i=0; i<BCv_transp.size(); i++){
				if(BCv_transp[i].label=="DIR")
					BCv_transp[i].value = 0;
			}
	
		asm_network_bc_transp(Fv, Avv, mimv, mf_Cv, mf_coefv, BCv_transp, beta_v, param.R() );
		}
//===========================================================================
	gmm::add(Avv, AM_temp_v);
	gmm::add(Fv, FM_temp_v);
	// De-allocate memory
	gmm::clear(Avv);
	gmm::clear(Fv);
	
	}/* end of assembly_rhs_transp */


	void transport3d1d::update_network (vector_type Pigreco){

		/*
		At each time step, the right hand side is modified by the time derivative term.
		Since we must ensure in a strong way the Dirichlet conditions, by modifying the monolithic matrix and the rhs vector, we save both AM_transp and FM_transp, where are assembled the stationary terms; 	then, we work on AM_temp and FM_temp, modifying them when necessary.
		*/

		#ifdef M3D1D_VERBOSE_
			cout << "  Update monolithic matrix and rhs vector ..." << endl;
		#endif

		gmm::copy(AM_transp_v, AM_temp_v);
		gmm::copy(FM_transp_v, FM_temp_v);

		bool SATURATION  = PARAM.int_value("SATURATION","Flag to use the saturation per the adhesive term");
		
		if(SATURATION==1){
			sparse_matrix_type Adhv(dof_transp.Cv(), dof_transp.Cv()); gmm::clear(Adhv);
			asm_network_nano_transp(Adhv, mimv, mf_Cv, mf_coefv, Pigreco); 	//if Pigreco is of type mf_coefv		
			masslumping(Adhv);
			gmm::add(Adhv, AM_temp_v);
		}

		// update rhs (time step mass term)
		vector_type TFv(dof_transp.Cv());
		vector_type RR(dof.coefv()); gmm::clear(RR);
		gmm::add(param.R(),RR); gmm::vscale(param.R(),RR); gmm:: scale (RR, pi);
		asm_source_term(TFv,mimv, mf_Cv, mf_coefv, RR);
		gmm::vscale(UM_transp_v, TFv);
		gmm::scale(TFv, (1.0/param_transp.dt())); // dt time step
		gmm::add(TFv, FM_temp_v);
		gmm::clear(TFv);

		assembly_rhs_network_transp();
	
	} /* end of update_transp*/

	
	void transport3d1d::update_tissue (vector_type RHS_nano){

		/*
		At each time step, the right hand side is modified by the time derivative term.
		Since we must ensure in a strong way the Dirichlet conditions, by modifying the monolithic matrix and the rhs vector, we save both AM_transp and FM_transp, where are assembled the stationary terms; 	then, we work on AM_temp and FM_temp, modifying them when necessary.
		*/

		#ifdef M3D1D_VERBOSE_
			cout << "  Update monolithic matrix and rhs vector ..." << endl;
		#endif

		gmm::copy(AM_transp_t, AM_temp_t);
		gmm::copy(FM_transp_t, FM_temp_t);

		// update rhs (time step mass term)
		vector_type TFt(dof_transp.Ct());
		asm_source_term(TFt,mimt, mf_Ct, mf_Ct,UM_transp_t); 
		gmm::scale(TFt, (1.0/param_transp.dt())); // dt time step
		gmm::add(TFt, FM_temp_t);
		gmm::clear(TFt); 

		gmm::scale(RHS_nano, 1/(param_transp.cNP()*param_transp.Uadim())); 	
		gmm::add(RHS_nano, FM_temp_t);	

		//update rhs (bundary condition terms)
		assembly_rhs_tissue_transp();

	} /* end of update_transp*/

	
	bool transport3d1d::solve_transp (void)
 	{
	  	#ifdef M3D1D_VERBOSE_
		cout << "Solving the monolithic system ... " << endl;
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


		vector_type Cvold(dof_transp.Cv());gmm::clear(Cvold);
		vector_type Cvnew(dof_transp.Cv()); gmm::clear(Cvnew);

		vector_type Psiold(mf_coefv.nb_dof()); gmm::clear(Psiold);
		gmm::resize(Psi, mf_coefv.nb_dof()); gmm::clear(Psi);
		
		vector_type psi_max_vec(mf_coefv.nb_dof()); gmm::clear(psi_max_vec);
		vector_type Pigrecostar(mf_coefv.nb_dof()); gmm::clear(Pigrecostar);
		 
		if(NANO==1 && SATURATION==1)
			psi_max_vec.assign(mf_coefv.nb_dof(),  param_transp.psi_max());
		
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

		const bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
		asm_exchange_aux_mat_transp(Mbar1, Mlin1, 
					mimv, mf_Ct, mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);

		// bluid oncotic term  ---> it is 0 for the drug delivery formulation
		vector_type ONCOTIC (dof.Pv()); gmm::clear(ONCOTIC);

		gmm::mult_add(gmm::scaled(Mbar1,-1.0), 
			  gmm::sub_vector(UM, gmm::sub_interval(dof.Ut(), dof.Pt())), ONCOTIC);

		scalar_type picoef=param.sigma()*(param.pi_v()-param.pi_t());
	    vector_type DeltaPi(dof.Pv(),picoef);
	    gmm::add(gmm::scaled(DeltaPi,-1.0), ONCOTIC);	
		gmm::scale(ONCOTIC,0.5*(1.0-param.sigma())*param.Q(0));

		// We set the oncotic term to 0 because we need Btv just with coefficient 2piR to build the source term
		gmm::scale(ONCOTIC,0);			
		// build permeability term  ---> it is just 2\pi R for drug delivery formulation
		vector_type PERM (dof.coefv());
		gmm::copy(param.R(), PERM);
		//gmm::scale(PERM, 2*pi*param_transp.Y()[0]);
		//We set thr PERM just equal to 2piR
		gmm::scale(PERM, 2*pi);

		//build exchange matrixes	
		asm_exchange_mat_transp(Btt1, Btv1, Bvt1, Bvv1,
				mimv, mf_Cv, mf_coefv, mf_Pv, Mbar1, Mlin1, 
				ONCOTIC, PERM, NEWFORM);

		std::ofstream integrated_ct_file(descr_transp.OUTPUT+"Ct_integrated.txt"); 

		for(t=0;t<=param_transp.T()*(!descr_transp.STATIONARY) ; t = t + param_transp.dt() + (param_transp.dt()==0)){
			time_count++; 
			scalar_type time_partial=gmm::uclock_sec();

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

				//=====================================================================================================
				//computing Psimean for each branch and at each time instant
				string time_suff_ = "";
				std::ostringstream convert_;
				convert_ << time_count;
				time_suff_ = convert_.str();

				for(size_type i=0; i<mf_coefv.nb_dof(); i++) { 
	      			Pigrecostar[i] = PiGrecoproj_scaled[i]*std::max(0.0,(psi_max_vec[i]-Psi[i]));
	      		}
				gmm::scale(Pigrecostar, 1.0/param_transp.psi_max());
			
				gmm::clear(Psiold);
				gmm::add(Psi, Psiold);
			
			}	

			//Print Psi at several time steps
			string time_suff = "";
			std::ostringstream convert;
			convert << time_count;
			time_suff = convert.str();
			
			cout << "  Exporting Psi ..." << endl;
			vtk_export exp_Psi(descr_transp.OUTPUT+"Psi_t"+time_suff+".vtk");
			exp_Psi.exporting(mf_coefv);
			exp_Psi.write_mesh();
			exp_Psi.write_point_data(mf_coefv, Psi, "Psi");

			
			if(descr_transp.STATIONARY){
				std::cout<<"Stationary problem... "<<std::endl;	
			}
			else{
				std::cout<<"-------------------------------------------"<<std::endl;
				std::cout<<"Iteration number: "<<time_count<<std::endl;
				std::cout<<"time = "<<t<<" s"<<std::endl<<std::endl;	
			}

			//Update rhs and boundary condition
			update_network(Pigrecostar);
			
			gmm::csc_matrix<scalar_type> A_transp_v;
			gmm::clean(AM_transp_v, 1E-12);
			gmm::copy(AM_temp_v, A_transp_v);
			
			vector_type F_transp_v(gmm::vect_size(FM_transp_v));
			gmm::clean(FM_transp_v, 1E-12);
			gmm::copy(FM_temp_v, F_transp_v);
			
			#ifdef M3D1D_VERBOSE_
				cout << "  Applying the SuperLU method to network problem... " << endl;
			#endif
			scalar_type cond_v;
			gmm::SuperLU_solve(A_transp_v, UM_transp_v, F_transp_v, cond_v);
			cout << "  Condition number (transport network problem): " << cond_v << endl;
			
					
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

			//Scalo Psi per la concentrazione di riferimento			
			gmm::copy(gmm::scaled(Psi_proj,2.39E14), J);

			scalar_type J_coeff;
			if(t==0)
				J_coeff=0.0;
			else{
				J_coeff= param_transp.ups()*param_transp.b()*pow(t,param_transp.b()-1.0)*param_transp.cNP()*param_transp.VNP()/pow((pow(t,param_transp.b())+param_transp.ups()),2.0);
			}
				
			gmm::scale(J, J_coeff);	
			vector_type RHS_nano_(dof_transp.Ct()); gmm::clear(RHS_nano_);
			gmm::mult(gmm::scaled(Btv1,-1.0),J,RHS_nano_);

			update_tissue(RHS_nano_);
				
			gmm::csc_matrix<scalar_type> A_transp_t;
			gmm::clean(AM_transp_t, 1E-12);
			gmm::copy(AM_temp_t, A_transp_t);
			
			vector_type F_transp_t(gmm::vect_size(FM_transp_t));
			gmm::clean(FM_transp_t, 1E-12);
			gmm::copy(FM_temp_t, F_transp_t);
			
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the SuperLU method to tissue problem... " << endl;
			#endif
			scalar_type cond_t;
			gmm::SuperLU_solve(A_transp_t, UM_transp_t, F_transp_t, cond_t);
			cout << "  Condition number (transport tissue problem): " << cond_t << endl;


			// computing the integral of Ct over time
			integrated_ct_file<< asm_integral(mf_Ct, mimt, UM_transp_t)<< endl;
			
			time_suff = "";
			std::ostringstream convert_;
			convert_ << time_count;
			time_suff = convert_.str();
			export_vtk_transp(time_suff); 


			#ifdef M3D1D_VERBOSE_		
				std::cout<<"exported!"<<std::endl;
			#endif	
			
			if(!descr_transp.STATIONARY)
				cout << "... time to solve : "	<< gmm::uclock_sec() - time_partial << " seconds\n";
		
		} //end of cycle over time 
		if(!descr_transp.STATIONARY){
			cout << endl<<"... time to solve all the time steps: " << gmm::uclock_sec() - time << " seconds\n";				}
		else{
			cout << endl<<"... time to solve : " << gmm::uclock_sec() - time << " seconds\n";
		}

		integrated_ct_file.close();
		
		//Add the last iteration to Psi
		if (NANO==1 && SATURATION==1){
			vector_type f_Cvold_proj(mf_coefv.nb_dof()); gmm::clear(f_Cvold_proj);
			vector_type Cvold_proj(mf_coefv.nb_dof()); gmm::clear(Cvold_proj);

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
			
			for(size_type i=0; i<mf_coefv.nb_dof(); i++) { 
					Psi[i] = PiGrecoproj[i]*std::max(0.0,(psi_max_vec[i]-Psiold[i]))*Cvold_proj[i];
			}
			gmm::scale(Psi, param_transp.dadim()*param_transp.dt()/param_transp.psi_max());
			gmm::add(Psiold, Psi);

			//Print Psi at several time steps
			string time_suff = "";
			std::ostringstream convert;
			convert << time_count;
			time_suff = convert.str();
			
			if(iteraz%5==0){ 
			cout << "  Exporting Psi ..." << endl;
			vtk_export exp_Psi(descr_transp.OUTPUT+"Psi_t"+time_suff+".vtk");
			exp_Psi.exporting(mf_coefv);
			exp_Psi.write_mesh();
			exp_Psi.write_point_data(mf_coefv, Psi, "Psi");
			}
		}
		

		export_vtk_nano();

		return true;
 	}; // end of solve_transp
	



	//Compute the residuals for mass balance at each junction 
	void transport3d1d::mass_balance(void){

	#ifdef M3D1D_VERBOSE_		
	cout << " Compute MBD and MBA "   << endl;
	#endif	

	// initialize the MBD and MBA to zero (clear at eac time step)
	for (size_type i=0; i<Jv_transp.size(); ++i){
		Jv_transp[i].MBD=0;
		Jv_transp[i].MBA=0;
	}	

	size_type shift = 0; //counter for branches
	
	for (size_type i=0; i<mf_Uvi.size(); ++i){ // branch loop 	
		if(i>0) shift += mf_Uvi[i-1].nb_dof(); 
		mesh_region &rg_branch = meshv.region(i); // branch region

		for (size_type j=0; j<Jv_transp.size(); ++j){ //junction loop
			mesh_region &rg_junction = meshv.region(Jv_transp[j].rg);  // junction region
			// Iterators for all the branches which flow in the junction j
			std::vector<long signed int>::const_iterator bb = Jv_transp[j].branches.begin();
			std::vector<long signed int>::const_iterator be = Jv_transp[j].branches.end();
		
			//Check if outflow of branch i is in junction j			
			if ((std::find(bb, be, +i) != be)){

				// find last element of the branch
				getfem::mr_visitor ii(rg_branch); int temp=0;
				for(; !ii.finished() && temp==0; ++ii){ // loop for convexes of the branch i
					getfem::mr_visitor ii_temp=ii;
					++ii_temp;
					if(ii_temp.finished()) temp=1;
				} 
				
				//import radius of the branch
				scalar_type Ri = param.R(mimv, i);
				// find the dof for last point for Cv e Uv
				size_type last_C, last_C2, last_U;
				vector_type dof_enum_C,dof_enum_U;
				int fine_C=0, fine_U=0;
				for (mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
				for (auto b : mf_Cv.ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_C.emplace_back(b);
					fine_C++;}			
				for (auto b : mf_Uvi[i].ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_U.emplace_back(b);
					fine_U++;}			
				}
	

				last_C = dof_enum_C[fine_C-1];
				last_C2= dof_enum_C[fine_C-2];
				last_U = dof_enum_U[fine_U-1];
				dof_enum_C.clear(); dof_enum_U.clear();
	
				scalar_type DIFF=0, ADV=0;

				//Compute the diffusive flux
				DIFF = pi* Ri * Ri*(
						UM_transp_v[last_C]-UM_transp_v[last_C2] )
						/estimate_h(meshv, ii.cv()) ;
				//Compute the advective fluxes
				ADV = pi* Ri * Ri*UM[dof.Ut()+dof.Pt()+shift+last_U]*UM_transp_v[last_C];
	
				#ifdef M3D1D_VERBOSE_		
				cout << "------------------------------------------ "   << endl;
				cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_transp[j].rg<<" )"<<endl;
				cout<<"MBD_partial = "<< DIFF<< endl;
				cout<<"MBA_partial = "<< ADV<<  endl;
				#endif	
		
				Jv_transp[j].MBD -= DIFF;
				Jv_transp[j].MBA -= ADV;
					
			}// end of check for outflow branch
	
			//Check if inflow of branch i is in junction j			
			if ( (i!=0) &&  (std::find(bb, be, -i) != be )){  //notice that, in build_vessel_transp, we assume that the branch zero cannot have the inflow in a junction. 


				// find first element of the branch
				getfem::mr_visitor ii(rg_branch);
			
				//import radius of the branch
				scalar_type Ri = param.R(mimv, i);
				// find the dof for last point for Cv e Uv
				size_type first_C, first_C2, first_U;
				vector_type dof_enum_C,dof_enum_U;
				int fine_C=0, fine_U=0;
				for (mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
				for (auto b : mf_Cv.ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_C.emplace_back(b);
					fine_C++;}			
				for (auto b : mf_Uvi[i].ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_U.emplace_back(b);
					fine_U++;}			
				}
	
				first_C = dof_enum_C[0];
				first_C2= dof_enum_C[1];
				first_U = dof_enum_U[0];
				dof_enum_C.clear(); dof_enum_U.clear();
	
				scalar_type DIFF=0,ADV=0;
	
				//Compute the diffusive flux
				 DIFF = pi* Ri * Ri*(
						UM_transp_v[first_C2]-UM_transp_v[first_C] )
						/estimate_h(meshv, ii.cv()) ;
	
	
				//Compute the advective fluxes
				ADV = pi* Ri * Ri*UM[dof.Ut()+dof.Pt()+shift+first_U]*UM_transp_v[first_C];
	
				#ifdef M3D1D_VERBOSE_		
				cout << "------------------------------------------ "   << endl;
				cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_transp[j].rg<<" )"<<endl;
				cout<<"MBD_partial = "<< DIFF<< endl;
				cout<<"MBA_partial = "<< ADV<<  endl;
				#endif	
		
				Jv_transp[j].MBD += DIFF;
				Jv_transp[j].MBA += ADV;
					
			}// end of check for outflow branch

		} //end of junction loop
	} // end of branch loop	





	cout << "  Junctions: " << endl;
	for (size_type i=0; i<Jv_transp.size(); ++i){
		cout << "    -  label=" << Jv_transp[i].label 
			 << ", value=" << Jv_transp[i].value << ", ind=" << Jv_transp[i].idx 
			 << ", rg=" << Jv_transp[i].rg << ", branches=" << Jv_transp[i].branches << endl; 
		cout << " Mass balance of diffusive fluxes = " << Jv_transp[i].MBD << endl; 
		cout << " Mass balance of advective fluxes = " << Jv_transp[i].MBA << endl;
		cout << "             ------------------- "   << endl;
	} 	
		cout << "----------------------------------------------- "   << endl;

	
	}; // end of mass_balance


 

 const void transport3d1d::export_vtk_transp (const string & time_suff,const string & suff)
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
	cout << "  Exporting Ct ..." << endl;
	#endif
	vtk_export exp_Ct(descr_transp.OUTPUT+"Ct"+suff+"_t"+time_suff+".vtk");
	exp_Ct.exporting(mf_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct, Ct, "Ct");

	// Writing on txt file
	std::ofstream Ct_file(descr_transp.OUTPUT+"Ct_c"+suff+"_t"+time_suff+".txt"); 
	Ct_file<<"BEGIN_LIST" << std::endl;
	for(scalar_type i=0;i<mf_Ct.nb_dof();i++){
		Ct_file<<Ct[i];
		Ct_file<<std::endl;
	}
	Ct_file<< "END_LIST" << std::endl;



	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+"Cv"+suff+"_t"+time_suff+".vtk");
	exp_Cv.exporting(mf_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_Cv, Cv, "Cv");

	// Writing on txt file
	std::ofstream Cv_file(descr_transp.OUTPUT+"Cv_c"+suff+"_t"+time_suff+".txt"); 
	Cv_file<<"BEGIN_LIST" << std::endl;
	for(scalar_type i=0;i<mf_Cv.nb_dof();i++){
		Cv_file<<Cv[i];
		Cv_file<<std::endl;
	}
	Cv_file<< "END_LIST" << std::endl;

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif

//	If you don't want to use vtk:
/*	use this method to export a vector:

	std::ofstream outFF(descr_transp.OUTPUT+"FF.txt");
		outFF << gmm::col_vector(F_transp);
		outFF.close(); 
*/ 

/*	use this method to export a matrix 

	gmm::MatrixMarket_IO::write("SM.mm" , SM);
*/

  }
 }; // end of export_transp
 
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

	// Writing on txt file
	std::ofstream Re_file(descr_transp.OUTPUT+"Re_c.txt"); 
	Re_file<<"BEGIN_LIST" << std::endl;
	for(scalar_type i=0;i<mf_coefv.nb_dof();i++){
		Re_file<<Reproj[i];
		Re_file<<std::endl;
	}
	Re_file<< "END_LIST" << std::endl;

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Probability of adhesion_v..." << endl;
	#endif

	vtk_export exp_padhv(descr_transp.OUTPUT+"Padhproj.vtk");
	exp_padhv.exporting(mf_coefv);
	exp_padhv.write_mesh();
	exp_padhv.write_point_data(mf_coefv, P_adhproj, "Padhv"); 
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Wall shear rate_v ..." << endl;
	#endif
		
	vtk_export exp_WSRv(descr_transp.OUTPUT+"WSRproj.vtk");
	exp_WSRv.exporting(mf_coefv);
	exp_WSRv.write_mesh();
	exp_WSRv.write_point_data(mf_coefv, WSRproj, "WSRv"); 

	// Writing on txt file
	std::ofstream WSR_file(descr_transp.OUTPUT+"WSR_c.txt"); 
	WSR_file<<"BEGIN_LIST" << std::endl;
	for(scalar_type i=0;i<mf_coefv.nb_dof();i++){
		WSR_file<<WSRproj[i];
		WSR_file<<std::endl;
	}
	Re_file<< "END_LIST" << std::endl;
		
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting vascular adhesion parameter Pi ..." << endl;
	#endif
		
	vtk_export exp_pigrecov(descr_transp.OUTPUT+"PiGrecoproj.vtk");
	exp_pigrecov.exporting(mf_coefv);
	exp_pigrecov.write_mesh();
	exp_pigrecov.write_point_data(mf_coefv, PiGrecoproj, "PiGreco_adim"); 

		#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Uv abs ..." << endl;
	#endif
		
	vtk_export exp_uv(descr_transp.OUTPUT+"Uvproj_abs.vtk");
	exp_uv.exporting(mf_coefv);
	exp_uv.write_mesh();
	exp_uv.write_point_data(mf_coefv, Uv_abs, "Uvproj_abs"); 

	//computing the total mean velocity 
	scalar_type Uv_mean = asm_integral(mf_coefv, mimv, Uv_abs);
	getfem::generic_assembly assem1;
        assem1.set("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
        assem1.push_mi(mimv);
        assem1.push_mf(mf_coefv);
        std::vector<scalar_type> one(gmm::vect_size(Uv_abs), 1.0);
        assem1.push_data(one);
        std::vector<scalar_type> measure(1);
        assem1.push_vec(measure);
        assem1.assembly();
	
     Uv_mean = Uv_mean/measure[0];

     //computing Psi mean
     scalar_type Psi_network_mean = asm_integral(mf_coefv, mimv, Psi);
     std::cout << "Psi integrated " << Psi_network_mean << std::endl;
     vtk_export exp_(descr_transp.OUTPUT+"prova.vtk");
	exp_.exporting(mf_coefv);
	exp_.write_mesh();
	exp_.write_point_data(mf_coefv, Psi, "prova"); 
	std::cout << "measure[0] " << measure[0] << std::endl;
     Psi_network_mean = Psi_network_mean /measure[0];



		#ifdef M3D1D_VERBOSE_
	cout << "  Exporting WSS ..." << endl;
	#endif
	vector_type WSS(mf_coefv.nb_dof()); gmm::clear(WSS);
	for (size_type k=0; k<mf_coefv.nb_dof(); k++)
		WSS[k] = WSRproj[k]* MU[k];
	vtk_export exp_wss(descr_transp.OUTPUT+"WSS.vtk");
	exp_wss.exporting(mf_coefv);
	exp_wss.write_mesh();
	exp_wss.write_point_data(mf_coefv, WSS, "WSS"); 

	// Writing on txt file
	std::ofstream WSS_file(descr_transp.OUTPUT+"WSS_c.txt"); 
	WSS_file<<"BEGIN_LIST" << std::endl;
	for(scalar_type i=0;i<mf_coefv.nb_dof();i++){
		WSS_file<<WSS[i];
		WSS_file<<std::endl;
	}
	WSS_file<< "END_LIST" << std::endl;

	 /*    //summing up the dimensional variables characterizing the problem
     std::cout << "------------------------------------------------------------------------" << std::endl; 
     std::cout << "		diameter = [ " 
     	<< *std::min_element(param.R().begin(), param.R().end()) * 2 * param_transp.dadim() 
     	<< " , "
     	<<*std::max_element(param.R().begin(), param.R().end()) * 2 * param_transp.dadim()
     	<< " ]"
     	<< std::endl;

     std::cout << "		\nvelocity = [ " 
     	<< *std::min_element(Uv_abs.begin(), Uv_abs.end()) * param_transp.Uadim() *1E3
     	<< " , "
     	<<*std::max_element(Uv_abs.begin(), Uv_abs.end())  * param_transp.Uadim() *1E3
     	<< " ] mm/s"
     	<< std::endl;

     std::cout << "		\nmean velocity = " << Uv_mean * param_transp.Uadim() *1E3 << " mm/s"<< std::endl;

     std::cout << "		\nreynolds = [ " 
     	<< *std::min_element(Reproj.begin(), Reproj.end()) 
     	<< " , "
     	<<*std::max_element(Reproj.begin(), Reproj.end()) 
     	<< " ]"
     	<< std::endl;   

     std::cout << "		\nWSS = [ " 
     	<< *std::min_element(WSS.begin(), WSS.end()) * 10 
     	<< " , "
     	<<*std::max_element(WSS.begin(), WSS.end()) * 10
     	<< " ] dyn/cm^2"
     	<< std::endl;   

     std::cout << "		\nviscosity = [ " 
     	<< *std::min_element(MU.begin(), MU.end()) *1E3
     	<< " , "
     	<<*std::max_element(MU.begin(), MU.end()) *1E3
     	<< " ] cPS"
     	<< std::endl;  

     std::cout << "------------------------------------------------------------------------" << std::endl; 
*/
     std::ofstream resume(descr_transp.OUTPUT+"results.txt"); 
         //summing up the dimensional variables characterizing the problem
     resume << "------------------------------------------------------------------------" << std::endl; 
     resume << "		diameter = [ " 
     	<< *std::min_element(param.R().begin(), param.R().end()) * 2 * param_transp.dadim() 
     	<< " , "
     	<<*std::max_element(param.R().begin(), param.R().end()) * 2 * param_transp.dadim()
     	<< " ]"
     	<< std::endl;

     resume << "		\nvelocity = [ " 
     	<< *std::min_element(Uv_abs.begin(), Uv_abs.end()) * param_transp.Uadim() *1E3
     	<< " , "
     	<<*std::max_element(Uv_abs.begin(), Uv_abs.end())  * param_transp.Uadim() *1E3
     	<< " ] mm/s"
     	<< std::endl;

     resume << "		\nmean velocity = " << Uv_mean * param_transp.Uadim() *1E3 << " mm/s"<< std::endl;

     resume << "		\nreynolds = [ " 
     	<< *std::min_element(Reproj.begin(), Reproj.end()) 
     	<< " , "
     	<<*std::max_element(Reproj.begin(), Reproj.end()) 
     	<< " ]"
     	<< std::endl;   

     resume<< "		\nWSS = [ " 
     	<< *std::min_element(WSS.begin(), WSS.end()) * 10 
     	<< " , "
     	<<*std::max_element(WSS.begin(), WSS.end()) * 10
     	<< " ] dyn/cm^2"
     	<< std::endl;   

     resume << "		\nviscosity = [ " 
     	<< *std::min_element(MU.begin(), MU.end()) *1E3
     	<< " , "
     	<<*std::max_element(MU.begin(), MU.end()) *1E3
     	<< " ] cPS"
     	<< std::endl; 
    	
    //dimensionalization: Psi * cref #/m^3 = Psi #/m^2 ---> * 1E-6 = Psi #/mm^2;
     resume << "		\nmean psi = " << Psi_network_mean * 2.39*1E8  << " #/mm^2"<< std::endl;


     std::cout << "------------------------------------------------------------------------" << std::endl; 


     resume.close();


	/*#ifdef M3D1D_VERBOSE_
	cout << "  Exporting density of nanoparticles adhering the wall Psi ..." << endl;
	#endif
		
	vtk_export exp_psi(descr_transp.OUTPUT+"Psi.vtk");
	exp_psi.exporting(mf_coefv);
	exp_psi.write_mesh();
	exp_psi.write_point_data(mf_coefv, Psi, "Psi"); 

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif*/
  }
 };
  

  // Interface with problem3d1d class
  	//! Initialize the problem
	void transport3d1d::init_fluid (int argc, char *argv[])
	{ 
	problem3d1d::init(argc, argv);
	};
	
	//! Assemble the problem
	void transport3d1d::assembly_fluid (void)
	{ 
	problem3d1d::assembly();
	};
	//! Solve the problem
	bool transport3d1d::solve_fluid (void)
	{ 
	return problem3d1d::solve();
	};	
	//! Export the solution
	const void transport3d1d::export_vtk_fluid (const string & suff)
	{ 
	problem3d1d::export_vtk(suff);
	};
	
	

void transport3d1d::mesh_test(void){

	// build coarser meshes

	std::cout << "Build coarser meshes for tissue and vessels" << std::endl;

	mesh mesht_c;
	mesh meshv_c;


	mesht_c.clear();

	cout << "Building the mesh for the tissue  "   << endl;
	string stc("GT='GT_PK(3,1)'; " 
				   "NSUBDIV = [11,8,6];"   
				   "ORG = [0,0,0];"  
				   "SIZES = [7.4,5,4]; "   
				   "NOISED = 0" ); 

	regular_mesh(mesht_c, stc);

	cout << "Importing the 1D mesh for the vessel COARSE"   << endl;
	std::ifstream ifs_c("../meshes/ratt93b_coarse/rattm93b_norm50.pts");
	GMM_ASSERT1(ifs_c.good(), "impossible to read from file " << "../meshes/ratt93b_coarse/rattm93b_norm50.pts");
	meshv_c.clear();

	vector<node> BCv_c;
	vector_size_type nb_v_c;

	import_pts_file(ifs_c, meshv_c, BCv_c, nb_v_c, "GT_PK(1,1)");

	ifs_c.close();

	// finer meshes
	std::cout << "Build finer meshes for tissue and vessels" << std::endl;

	mesh mesht_f;
	mesh meshv_f;

	mesht_f.clear();
	cout << "Building the mesh for the tissue  "   << endl;
	string stf("GT='GT_PK(3,1)'; " 
				   "NSUBDIV= [22,15,12];" 
				   "ORG= [0,0,0];"   
				   "SIZES= [7.4,5,4] ; " 
				   "NOISED= 0" ); 
	regular_mesh(mesht_f, stf);

	cout << "Importing the 1D mesh for the vessel FINE"   << endl;
	std::ifstream ifs_f("../meshes/ratt93b/rattm93b_norm50.pts");
	GMM_ASSERT1(ifs_f.good(), "impossible to read from file " << "../meshes/ratt93b/rattm93b_norm50.pts");

	vector<node> BCv_f;
	vector_size_type nb_v_f;	
	meshv_f.clear();
	import_pts_file(ifs_f, meshv_f, BCv_f, nb_v_f, "GT_PK(1,1)");

	ifs_f.close();


// Building the mesh fem

	std::cout << "Setting IMs for tissue and vessel problems on the COARSE mesh..." << std::endl;

	pintegration_method pim_t = int_method_descriptor(descr.IM_TYPET);
	pintegration_method pim_v = int_method_descriptor(descr.IM_TYPEV);

	bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(descr.MESH_TYPET);
	bgeot::pgeometric_trans pgt_v = bgeot::geometric_trans_descriptor(descr.MESH_TYPEV);

	pfem pf_Pt = fem_descriptor(descr.FEM_TYPET_P);
	pfem pf_Pv = fem_descriptor(descr.FEM_TYPEV_P);

	//DIMT = pgt_t->dim();	//DIMV = 1;

	mesh_im mimt_c(mesht_c), mimv_c(meshv_c);
	mimt_c.set_integration_method(mesht_c.convex_index(), pim_t);
	mimv_c.set_integration_method(meshv_c.convex_index(), pim_v);
		
	mesh_fem mf_Pt_c(mesht_c), mf_Pv_c(meshv_c);
	//mf_Pt_c.set_qdim(bgeot::dim_type(DIMT)); 
	mf_Pt_c.set_finite_element(mesht_c.convex_index(), pf_Pt);
	mf_Pv_c.set_finite_element(meshv_c.convex_index(), pf_Pv);

	std::cout << "Setting IMs for tissue and vessel problems on the FINE mesh..." << std::endl;

	mesh_im mimt_f(mesht_f), mimv_f(meshv_f);
	mimt_f.set_integration_method(mesht_f.convex_index(), pim_t);
	mimv_f.set_integration_method(meshv_f.convex_index(), pim_v);
	
	mesh_fem mf_Pt_f(mesht_f), mf_Pv_f(meshv_f);
	//mf_Pt_f.set_qdim(bgeot::dim_type(DIMT)); 
	mf_Pt_f.set_finite_element(mesht_f.convex_index(), pf_Pt);
	mf_Pv_f.set_finite_element(meshv_f.convex_index(), pf_Pv);

	std::cout << "--------dof Pt_c " << mf_Pt_c.nb_dof() <<"--------dof Pt_f " <<  mf_Pt_f.nb_dof()<< std::endl;
	std::cout << "--------dof Pv_c " << mf_Pv_c.nb_dof() <<"--------dof Pv_f " <<  mf_Pv_f.nb_dof()<< std::endl;	

// Read the coarse solutions Pt and Pv -> save them in Pt_c Pv_c

	// Reading Pt on COARSE MESH
	vector_type Pt_c;
	std::ifstream ist("./coarse_txt/Pt_c.txt");
	if (!ist) cerr << "impossible to read from file " << "./coarse_txt/Pt_c.txt" << endl;
	ist.precision(16);
	ist.seekg(0); ist.clear();
	GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"), "This seems not to be a data file");
	std::string line;
	bool thend = false; 
	while (!thend){
		bgeot::get_token(ist, line, 1023);
		thend = (line=="END_LIST");
		if (!thend){
			Pt_c.emplace_back(stof(line));
		} 
	}	

	std::cout << "length of Pt_c " << Pt_c.size()<< " dof mf_Pt_c " << mf_Pt_c.nb_dof() << std::endl; 

	vtk_export exp_Pt_c(descr.OUTPUT+"Pt_c.vtk");
	exp_Pt_c.exporting(mf_Pt_c);
	exp_Pt_c.write_mesh();
	exp_Pt_c.write_point_data(mf_Pt_c, Pt_c, "Pt_c");


	// Reading Pv on COARSE MESH
	vector_type Pv_c;
	std::ifstream ist1("./coarse_txt/Pv_c.txt");
	if (!ist1) cerr << "impossible to read from file " << "./coarse_txt/Pv_c.txt" << endl;
	ist1.precision(16);
	ist1.seekg(0); ist1.clear();

	GMM_ASSERT1(bgeot::read_until(ist1, "BEGIN_LIST"), "This seems not to be a data file");
	thend = false; 
	while (!thend){
		bgeot::get_token(ist1, line, 1023);
		thend = (line=="END_LIST");
		if (!thend){
			Pv_c.emplace_back(stof(line));
		} 
	}	

	std::cout << "length of Pv_c " << Pv_c.size()<< " dof mf_Pv_c " << mf_Pv_c.nb_dof() << std::endl; 

	vtk_export exp_Pv_c(descr.OUTPUT+"Pv_c.vtk");
	exp_Pv_c.exporting(mf_Pv_c);
	exp_Pv_c.write_mesh();
	exp_Pv_c.write_point_data(mf_Pv_c, Pv_c, "Pv_c");


// Read Pt and Pv on the FINE mesh -> save them in Pv_f Pt_f

		// Reading Pt on FINE MESH
	vector_type Pt_f;
	std::ifstream ist2("./fine_txt/Pt_c.txt");
	if (!ist2) cerr << "impossible to read from file " << "./fine_txt/Pt_c.txt" << endl;
	ist2.precision(16);
	ist2.seekg(0); ist2.clear();
	GMM_ASSERT1(bgeot::read_until(ist2, "BEGIN_LIST"), "This seems not to be a data file");
	thend = false; 
	while (!thend){
		bgeot::get_token(ist2, line, 1023);
		thend = (line=="END_LIST");
		if (!thend){
			Pt_f.emplace_back(stof(line));
		} 
	}	

	std::cout << "length of Pt_f " << Pt_f.size()<< " dof mf_Pt_f " << mf_Pt_f.nb_dof() << std::endl; 

	vtk_export exp_Pt_f(descr.OUTPUT+"Pt_f.vtk");
	exp_Pt_f.exporting(mf_Pt_f);
	exp_Pt_f.write_mesh();
	exp_Pt_f.write_point_data(mf_Pt_f, Pt_f, "Pt_f");


	// Reading Pv on FINE MESH
	vector_type Pv_f;
	std::ifstream ist3("./fine_txt/Pv_c.txt");
	if (!ist3) cerr << "impossible to read from file " << "./fine_txt/Pv_c.txt" << endl;
	ist3.precision(16);
	ist3.seekg(0); ist3.clear();
	GMM_ASSERT1(bgeot::read_until(ist3, "BEGIN_LIST"), "This seems not to be a data file");
	thend = false; 
	while (!thend){
		bgeot::get_token(ist3, line, 1023);
		thend = (line=="END_LIST");
		if (!thend){
			Pv_f.emplace_back(stof(line));
		} 
	}	

	std::cout << "length of Pv_f " << Pv_f.size()<< " dof mf_Pv_f " << mf_Pv_f.nb_dof() << std::endl; 

	vtk_export exp_Pv_f(descr.OUTPUT+"Pv_f.vtk");
	exp_Pv_f.exporting(mf_Pv_f);
	exp_Pv_f.write_mesh();
	exp_Pv_f.write_point_data(mf_Pv_f, Pv_f, "Pv_f");


// Projection of Pt_c on Pt_f

	//creating the interpolated mesh fem on the fine mesh
	getfem::mesh_fem mf_interpole(mf_Pt_f.linked_mesh());
	pfem ifem = getfem::new_interpolated_fem(mf_Pt_c, mimt_f);
	dal::bit_vector nn = mf_Pt_f.convex_index();
	mf_interpole.set_finite_element(nn, ifem);

	/*//assembling the right hand side for projection
	vector_type rhs_Pt_proj(mf_Pt_f.nb_dof()); gmm::clear(rhs_Pt_proj);
	getfem::generic_assembly
	assempt("Pc=data$1(#1);" 
		"V$1(#2)+=comp(Base(#1).Base(#2))(j,:).Pc(j);"
		);
	assempt.push_mi(mimt_f);
	assempt.push_mf(mf_interpole);
	assempt.push_mf(mf_Pt_f);
	assempt.push_data(Pt_c);
	assempt.push_vec(rhs_Pt_proj);
	assempt.assembly();   //region 0 is the face x=0	
	*/
	// Computing the residual in the L2 norm

	std::vector<scalar_type> L2v(1);

	getfem::generic_assembly
	assemL2("Pc=data$1(#2);" "Pf=data$2(#1);"
		"t1=comp(Base(#2).Base(#2))(i,j).Pc(i).Pc(j);"
		"t2=comp(Base(#1).Base(#1))(i,j).Pf(i).Pf(j);"
		"t3=comp(Base(#2).Base(#1))(i,j).Pc(i).Pf(j);"
	"V$1() += t1 + t2 - t3-t3;");
	assemL2.push_mi(mimt_f);
	assemL2.push_mf(mf_Pt_f);
	assemL2.push_mf(mf_interpole);
	assemL2.push_data(Pt_c);
	assemL2.push_data(Pt_f);
	assemL2.push_vec(L2v);
	assemL2.assembly();   //region 0 is the face x=0	

	std::cout << "-----normL2(Pt_c - Pt_f) = " << sqrt(L2v[0]) << std::endl;

	vtk_export exp_Pproj(descr.OUTPUT+"Pt_proj.vtk");
	exp_Pproj.exporting(mf_interpole);
	exp_Pproj.write_mesh();
	exp_Pproj.write_point_data(mf_interpole, Pt_c, "Pt_c");

	del_interpolated_fem(ifem);

	/*
	//mass matrix for projection
	sparse_matrix_type M_proj(mf_Pt_f.nb_dof(), mf_Pt_f.nb_dof()); gmm::clear(M_proj);
	asm_mass_matrix(M_proj, mimt_f, mf_Pt_f );

	// solve the projection
	scalar_type cond;
	//Reynolds
	vector_type Pt_proj(mf_Pt_f.nb_dof());
	gmm::SuperLU_solve(M_proj, Pt_proj, rhs_Pt_proj,cond);
	
	vtk_export exp_Pproj(descr.OUTPUT+"Pt_proj.vtk");
	exp_Pproj.exporting(mf_Pt_f);
	exp_Pproj.write_mesh();
	exp_Pproj.write_point_data(mf_Pt_f, Pt_proj, "Pt_f");


	vector_type diff (mf_Pt_f.nb_dof());
	gmm::add( Pt_proj, diff);
	gmm::add(diff, gmm::scaled(Pt_f,-1.0), diff);
	scalar_type norm=0;
	scalar_type norm_Pt_f=0;
	for (size_type i =0 ; i< diff.size(); i++){
		norm += diff[i]*diff[i];
		norm_Pt_f +=  Pt_f[i]*Pt_f[i];
	}
	std::cout << " norm L2 vectors Pt_c and Pt_f" << sqrt(norm) << std::endl;
	std::cout << " norm  Pt_f" << sqrt(norm_Pt_f) << std::endl;
	std::cout << " realtive norm L2 vectors Pt_c and Pt_f" << sqrt(norm)/sqrt(norm_Pt_f) << std::endl;
	*/


} //end mesh_test
	
 
 } // end of namespace