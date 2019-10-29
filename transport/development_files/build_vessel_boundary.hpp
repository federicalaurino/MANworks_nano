void 
transport3d1d::build_vessel_boundary_transp2(void)
{ 
	#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
	#endif
try {

	dal::bit_vector extrema;   // global idx of extreme vertices in meshv
	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector junctions_region; // global idx of junctions regions in meshv

	Jv_transp.clear();
	nb_extrema=0; 
	nb_junctions=0;
	
	size_type fer = nb_branches; // First Empty Region
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
		if (meshv.convex_to_point(i0).size()==1){ 				/* inflow extremum */
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
		else if (meshv.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
			// Check if junction has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i0);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i0);
				junctions_region.add(fer);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 1); // single-face region
				// Create a new junction node
				Jv_transp.emplace_back("JUN", 0, i0, fer);
				fer++;
			}
			else{
			size_type fer_temp=0;
			while(junctions[fer_temp]!=i0)	fer_temp++;
			meshv.region(junctions_region[fer_temp]).add(cv,1);

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
					junctions_region.add(i1);
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
			else{
				size_type fer_temp=0;
				while(junctions[fer_temp]!=i1)	fer_temp++;
				meshv.region(junctions_region[fer_temp]).add(cv,0);}
			}
		}
		else if (meshv.convex_to_point(i1).size()>2){ /* non-trivial outflow junction */

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
				junctions_region.add(i1);
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

			
				size_type fer_temp=0;
				while(junctions[fer_temp]!=i1)	fer_temp++;
				meshv.region(junctions_region[fer_temp]).add(cv,0);

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


void transport3d1d::compute_mass_balance2(void){




	for (size_type i=0; i<Jv_transp.size(); ++i){
		Jv_transp[i].MBD=0;
		Jv_transp[i].MBA=0;
}
	size_type shift = 0;
	
	for (size_type i=0; i<mf_Uvi.size(); ++i){ // branch loop 
		
		if(i>0) shift += mf_Uvi[i-1].nb_dof(); 
	
		for (size_type j=0; j<Jv_transp.size(); ++j){ //junction loop

			// Iterators for all the branches which flow in the junction j
			std::vector<long signed int>::const_iterator bb = Jv_transp[j].branches.begin();
			std::vector<long signed int>::const_iterator be = Jv_transp[j].branches.end();

//Check if outflow of branch i is in junction j			
			if (std::find(bb, be, i) != be){

				//Import concentration value in network
				vector_type Cv(dof_transp.Cv()); 
				gmm::copy(gmm::sub_vector(UM_transp, 
					  gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);
				//import radius of te branch
				scalar_type Ri = param.R(mimv, i);
				// import diffusivity in te branch
				vector_type DIFF_param( dof.coefv());
				gmm::copy(param_transp.Av(), DIFF_param);
				gmm::scale(DIFF_param, pi*Ri*Ri);	
				// import velocity in the branch
 				vector_type ADV_param( mf_Uvi[i].nb_dof()); gmm::clear(ADV_param);
				gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  ADV_param);
				gmm::scale(ADV_param, pi*Ri*Ri);
				//Build auxiliary vectors in which we will temporarily store the fluxes
				//(generic assembly always needs vectors or matrixes, not scalars)
				std::vector<scalar_type> Diff(1);
				std::vector<scalar_type> Adv(1);

cout << "------------------------------------------ "   << endl;
cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_transp[j].rg<<" )"<<endl;

	generic_assembly 
	assem1("c=data$1(#1); l1=data$2(#2); l2=data$3(#2); l3=data$4(#2);  Diff=data$5(#3);"
		  "t=comp(Base(#3).Grad(#1).Base(#2));"
		  "V$1()+=t(d,i,1,j).Diff(d).c(i).l1(j)+t(d,i,2,j).Diff(d).c(i).l2(j)+t(d,i,3,j).Diff(d).c(i).l3(j);"); 
	assem1.push_mi(mimv);
	assem1.push_mf(mf_Cv);
	assem1.push_mf(mf_coefvi[i]);
	assem1.push_mf(mf_coefv);
	assem1.push_data(Cv);
	assem1.push_data(param.lambdax(i));
	assem1.push_data(param.lambday(i));
	assem1.push_data(param.lambdaz(i));
	assem1.push_data(DIFF_param);
	assem1.push_vec(Diff);
	assem1.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem2("c=data$1(#1); Adv=data$2(#2);"
		  "V()+= comp (Base(#1).Base(#2))(i,j).c(i).Adv(j);"); 
	assem2.push_mi(mimv);
	assem2.push_mf(mf_Cv);
	assem2.push_mf(mf_Uvi[i]);
	assem2.push_data(Cv);
	assem2.push_data(ADV_param);
	assem2.push_vec(Adv);
	assem2.assembly(Jv_transp[j].rg);

				Jv_transp[j].MBD -= Diff[0];
				Jv_transp[j].MBA -= Adv[0];

cout<<"MBD_partial = "<< Diff[0]<< "; MBD_cum = "<< Jv_transp[j].MBD <<";   DIFF_param = "<<DIFF_param[0]<<endl;
cout<<"MBA_partial = "<< Adv[0]<< "; MBA_cum = "<< Jv_transp[j].MBA <<";   ADV_param = "<<ADV_param[0]<<endl;

		}

			//Check if inflow of branch i is in junction j			
			if(i!=0)
			if ( std::find(bb, be, -i) != be){  //notice that, in build_vessel_transp, we assume that the branch zero cannot have the inflow in a junction. 


				//Import concentration value in network
				vector_type Cv(dof_transp.Cv()); 
				gmm::copy(gmm::sub_vector(UM_transp, 
					  gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);
				//import radius of te branch
				scalar_type Ri = param.R(mimv, i);
				// import diffusivity in te branch
				vector_type DIFF_param( dof.coefv());
				gmm::copy(param_transp.Av(), DIFF_param);
				gmm::scale(DIFF_param, pi*Ri*Ri);	
				// import velocity in the branch
 				vector_type ADV_param( mf_Uvi[i].nb_dof()); gmm::clear(ADV_param);
				gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  ADV_param);
				gmm::scale(ADV_param, pi*Ri*Ri);
				//Build auxiliary vectors in which we will temporarily store the fluxes
				//(generic assembly always needs vectors or matrixes, not scalars)
				std::vector<scalar_type> Diff(1);
				std::vector<scalar_type> Adv(1);

cout << "------------------------------------------ "   << endl;
cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_transp[j].rg<<" )"<<endl;

	generic_assembly 
	assem1("c=data$1(#1); l1=data$2(#2); l2=data$3(#2); l3=data$4(#2);  Diff=data$5(#3);"
		  "t=comp(Base(#3).Grad(#1).Base(#2));"
		  "V$1()+=t(d,i,1,j).Diff(d).c(i).l1(j)+t(d,i,2,j).Diff(d).c(i).l2(j)+t(d,i,3,j).Diff(d).c(i).l3(j);"); 
	assem1.push_mi(mimv);
	assem1.push_mf(mf_Cv);
	assem1.push_mf(mf_coefvi[i]);
	assem1.push_mf(mf_coefv);
	assem1.push_data(Cv);
	assem1.push_data(param.lambdax(i));
	assem1.push_data(param.lambday(i));
	assem1.push_data(param.lambdaz(i));
	assem1.push_data(DIFF_param);
	assem1.push_vec(Diff);
	assem1.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem2("c=data$1(#1); Adv=data$2(#2);"
		  "V()+= comp (Base(#1).Base(#2))(i,j).c(i).Adv(j);"); 
	assem2.push_mi(mimv);
	assem2.push_mf(mf_Cv);
	assem2.push_mf(mf_Uvi[i]);
	assem2.push_data(Cv);
	assem2.push_data(ADV_param);
	assem2.push_vec(Adv);
	assem2.assembly(Jv_transp[j].rg);

				Jv_transp[j].MBD -= Diff[0];
				Jv_transp[j].MBA -= Adv[0];


cout<<"MBD_partial = "<< Diff[0]<< "; MBD_cum = "<< Jv_transp[j].MBD <<";   DIFF_param = "<<DIFF_param[0]<<endl;
cout<<"MBA_partial = "<< Adv[0]<< "; MBA_cum = "<< Jv_transp[j].MBA <<";   ADV_param = "<<ADV_param[0]<<endl;


}

		} //end of junction loop
	} // end of branch loop	


	cout << "  Junctions: " << endl;
	for (size_type i=0; i<Jv_transp.size(); ++i){
		cout << "    -  label=" << Jv_transp[i].label 
			 << ", value=" << Jv_transp[i].value << ", ind=" << Jv_transp[i].idx 
			 << ", rg=" << Jv_transp[i].rg << ", branches=" << Jv_transp[i].branches << endl; 
		cout << " Mass balance of diffusive fluxes = " << Jv_transp[i].MBD << endl; 
		cout << " Mass balance of advective fluxes = " << Jv_transp[i].MBA << endl;
		cout << "			------------------- "   << endl;
	} 	cout << "---------------------------------------- "   << endl;
		



}; // end of compute_mass_balance



/*
	for (size_type i=0; i<Jv_transp.size(); ++i){
		Jv_transp[i].MBD=0;
		Jv_transp[i].MBA=0;
}
	size_type shift = 0;
	
	for (size_type i=0; i<mf_Uvi.size(); ++i){ // branch loop 
		
		if(i>0) shift += mf_Uvi[i-1].nb_dof(); 
	
		for (size_type j=0; j<Jv_transp.size(); ++j){ //junction loop

			// Iterators for all the branches which flow in the junction j
			std::vector<long signed int>::const_iterator bb = Jv_transp[j].branches.begin();
			std::vector<long signed int>::const_iterator be = Jv_transp[j].branches.end();
		
			//Check if outflow of branch i is in junction j			
			if (std::find(bb, be, i) != be){

				//Import concentration value in network
				vector_type Cv(dof_transp.Cv()); 
				gmm::copy(gmm::sub_vector(UM_transp, 
					  gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);
				//import radius of te branch
				scalar_type Ri = param.R(mimv, i);
				// import diffusivity in te branch
				vector_type DIFF_param( dof.coefv());
				gmm::copy(param_transp.Av(), DIFF_param);
				gmm::scale(DIFF_param, pi*Ri*Ri);	
				// import velocity in the branch
 				vector_type ADV_param( mf_Uvi[i].nb_dof()); gmm::clear(ADV_param);
				gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  ADV_param);
				gmm::scale(ADV_param, pi*Ri*Ri);
				//Build auxiliary vectors in which we will temporarily store the fluxes
				//(generic assembly always needs vectors or matrixes, not scalars)
				std::vector<scalar_type> DIFF(1);
				std::vector<scalar_type> ADV(1);

cout << "------------------------------------------ "   << endl;
cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_transp[j].rg<<" )"<<endl;

				asm_mass_fluxes(DIFF, ADV, 
					mimv, mf_Cv, mf_coefvi[i], mf_coefv, mf_Uvi[i],
					Cv, param.lambdax(i), param.lambday(i),param.lambdaz(i), 
					DIFF_param, ADV_param,					 
					mf_Cv.linked_mesh().region( Jv_transp[j].rg ) );

				Jv_transp[j].MBD -= DIFF[0];
				Jv_transp[j].MBA -= ADV[0];

std::vector<scalar_type> c(1);
std::vector<scalar_type> Diff(1);
std::vector<scalar_type> Adv(1);
std::vector<scalar_type> cv(3);
std::vector<scalar_type> u(1);
std::vector<scalar_type> a(1);
std::vector<scalar_type> l1(1);
std::vector<scalar_type> l2(1);
std::vector<scalar_type> l3(1);



	generic_assembly 
	assem3("c=data$1(#1);"
		  "V()+= comp (Base(#1))(i).c(i);"); 
	assem3.push_mi(mimv);
	assem3.push_mf(mf_Cv);
	assem3.push_data(Cv);
	assem3.push_vec(c);
	assem3.assembly(Jv_transp[j].rg);


	generic_assembly 
	assem4(" Adv=data$1(#1);"
		  "V()+= comp (Base(#1))(i).Adv(i);"); 
	assem4.push_mi(mimv);
	assem4.push_mf(mf_Uvi[i]);
	assem4.push_data(ADV_param);
	assem4.push_vec(u);
	assem4.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem5(" Diff=data$2(#2);"
		  "V()+= comp (Base(#2))(i).Diff(i);"); 
	assem5.push_mi(mimv);
	assem5.push_mf(mf_Cv);
	assem5.push_mf(mf_coefv);
	assem5.push_data(Cv);
	assem5.push_data(DIFF_param);
	assem5.push_vec(a);
	assem5.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem6("l1=data$2(#2); "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l1(j);"); 
	assem6.push_mi(mimv);
	assem6.push_mf(mf_Cv);
	assem6.push_mf(mf_coefvi[i]);
	assem6.push_mf(mf_coefv);
	assem6.push_data(Cv);
	assem6.push_data(param.lambdax(i));
	assem6.push_data(param.lambday(i));
	assem6.push_data(param.lambdaz(i));
	assem6.push_data(DIFF_param);
	assem6.push_vec(l1);
	assem6.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem7("l2=data$3(#2); "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l2(j);"); 
	assem7.push_mi(mimv);
	assem7.push_mf(mf_Cv);
	assem7.push_mf(mf_coefvi[i]);
	assem7.push_mf(mf_coefv);
	assem7.push_data(Cv);
	assem7.push_data(param.lambdax(i));
	assem7.push_data(param.lambday(i));
	assem7.push_data(param.lambdaz(i));
	assem7.push_data(DIFF_param);
	assem7.push_vec(l1);
	assem7.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem9("l3=data$4(#2);  "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l3(j);"); 
	assem9.push_mi(mimv);
	assem9.push_mf(mf_Cv);
	assem9.push_mf(mf_coefvi[i]);
	assem9.push_mf(mf_coefv);
	assem9.push_data(Cv);
	assem9.push_data(param.lambdax(i));
	assem9.push_data(param.lambday(i));
	assem9.push_data(param.lambdaz(i));
	assem9.push_data(DIFF_param);
	assem9.push_vec(l1);
	assem9.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem8("c=data$1(#1); "
		  "t=comp(Grad(#1));"
		  "V$1(3)+=t(i,:).c(i);"); 
	assem8.push_mi(mimv);
	assem8.push_mf(mf_Cv);
	assem8.push_mf(mf_coefvi[i]);
	assem8.push_mf(mf_coefv);
	assem8.push_data(Cv);
	assem8.push_data(param.lambdax(i));
	assem8.push_data(param.lambday(i));
	assem8.push_data(param.lambdaz(i));
	assem8.push_data(DIFF_param);
	assem8.push_vec(cv);
	assem8.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem1("c=data$1(#1); l1=data$2(#2); l2=data$3(#2); l3=data$4(#2);  Diff=data$5(#3);"
		  "t=comp(Base(#3).Grad(#1).Base(#2));"
		  "V$1()+=t(d,i,1,j).Diff(d).c(i).l1(j)+t(d,i,2,j).Diff(d).c(i).l2(j)+t(d,i,3,j).Diff(d).c(i).l3(j);"); 
	assem1.push_mi(mimv);
	assem1.push_mf(mf_Cv);
	assem1.push_mf(mf_coefvi[i]);
	assem1.push_mf(mf_coefv);
	assem1.push_data(Cv);
	assem1.push_data(param.lambdax(i));
	assem1.push_data(param.lambday(i));
	assem1.push_data(param.lambdaz(i));
	assem1.push_data(DIFF_param);
	assem1.push_vec(Diff);
	assem1.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem2("c=data$1(#1); Adv=data$2(#2);"
		  "V()+= comp (Base(#1).Base(#2))(i,j).c(i).Adv(j);"); 
	assem2.push_mi(mimv);
	assem2.push_mf(mf_Cv);
	assem2.push_mf(mf_Uvi[i]);
	assem2.push_data(Cv);
	assem2.push_data(ADV_param);
	assem2.push_vec(Adv);
	assem2.assembly(Jv_transp[j].rg);

cout<<"MBD partial calcolata in mass_balance = " << Diff[0] <<endl;
cout<<"MBA partial calcolata in mass_balance = " << Adv[0] <<endl;
cout<<"c = "<< c[0]<< "; dx(c) = "<< cv[0]<< "; dy(c) = "<< cv[1]<< "; dz(c) = "<< cv[2] <<"; "<<endl;
cout<<"diffusion = "<< a[0]<< "; advection = "<< u[0]<< "; l1 "<< l1[0]<< "; l2 = "<< l2[0]<< "; l3 = "<< l3[0] <<"; "<<endl;


vector_type Uv_basis( mf_Uvi[i].nb_dof()); gmm::clear(Uv_basis);
vector_type Cv_basis(dof_transp.Cv()); gmm::clear(Cv_basis);

generic_assembly
assem10( "V(#1)+=comp(Base(#1))(:);");
assem10.push_mi(mimv);
assem10.push_mf(mf_Cv);
assem10.push_vec(Cv_basis);
assem10.assembly(Jv_transp[j].rg);

generic_assembly
assem11("V(#1)+=comp(Base(#1))(:);");
assem11.push_mi(mimv);
assem11.push_mf(mf_Uvi[i]);
assem11.push_vec(Uv_basis);
assem11.assembly(Jv_transp[j].rg);


	string branch_suff = "";
	std::ostringstream convert;
	convert << i;
	branch_suff = convert.str();

	string junc_suff = "";
	std::ostringstream convert2;
	convert2 << j;
	junc_suff = convert2.str();

	std::ofstream outC("C_branch_"+branch_suff+"_junc_"+junc_suff+".txt");
		outC << gmm::col_vector(Cv_basis);
		outC.close();			

	std::ofstream outU("U_branch_"+branch_suff+"_junc_"+junc_suff+".txt");
		outU << gmm::col_vector(Uv_basis);
		outU.close();			


cout<<"MBD_partial = "<< DIFF[0]<< "; MBD_cum = "<< Jv_transp[j].MBD <<";   DIFF_param = "<<DIFF_param[0]<<endl;
cout<<"MBA_partial = "<< ADV[0]<< "; MBA_cum = "<< Jv_transp[j].MBA <<";   ADV_param = "<<ADV_param[0]<<endl;


				
}
			//Check if inflow of branch i is in junction j			
			if(i!=0)
			if ( std::find(bb, be, -i) != be){  //notice that, in build_vessel_transp, we assume that the branch zero cannot have the inflow in a junction. 


				//Import concentration value in network
				vector_type Cv(dof_transp.Cv()); 
				gmm::copy(gmm::sub_vector(UM_transp, 
					  gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);
				//import radius of te branch
				scalar_type Ri = param.R(mimv, i);
				// import diffusivity in te branch
				vector_type DIFF_param( dof.coefv());
				gmm::copy(param_transp.Av(), DIFF_param);
				gmm::scale(DIFF_param, pi*Ri*Ri);	
				// import velocity in the branch
 				vector_type ADV_param( mf_Uvi[i].nb_dof()); gmm::clear(ADV_param);
				gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  ADV_param);
				gmm::scale(ADV_param, pi*Ri*Ri);
				//Build auxiliary vectors in which we will temporarily store the fluxes
				//(generic assembly always needs vectors or matrixes, not scalars)
				std::vector<scalar_type> DIFF(1);
				std::vector<scalar_type> ADV(1);

cout << "------------------------------------------ "   << endl;
cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_transp[j].rg<<" )"<<endl;

				asm_mass_fluxes(DIFF, ADV, 
					mimv, mf_Cv, mf_coefvi[i], mf_coefv, mf_Uvi[i],
					Cv, param.lambdax(i), param.lambday(i),param.lambdaz(i), 
					DIFF_param, ADV_param,					 
					mf_Cv.linked_mesh().region( Jv_transp[j].rg ) );

				Jv_transp[j].MBD -= DIFF[0];
				Jv_transp[j].MBA -= ADV[0];


std::vector<scalar_type> c(1);
std::vector<scalar_type> Diff(1);
std::vector<scalar_type> Adv(1);
std::vector<scalar_type> cv(3);
std::vector<scalar_type> u(1);
std::vector<scalar_type> a(1);
std::vector<scalar_type> l1(1);
std::vector<scalar_type> l2(1);
std::vector<scalar_type> l3(1);



	generic_assembly 
	assem3("c=data$1(#1);"
		  "V()+= comp (Base(#1))(i).c(i);"); 
	assem3.push_mi(mimv);
	assem3.push_mf(mf_Cv);
	assem3.push_data(Cv);
	assem3.push_vec(c);
	assem3.assembly(Jv_transp[j].rg);


	generic_assembly 
	assem4(" Adv=data$2(#2);"
		  "V()+= comp (Base(#2))(i).Adv(i);"); 
	assem4.push_mi(mimv);
	assem4.push_mf(mf_Cv);
	assem4.push_mf(mf_Uvi[i]);
	assem4.push_data(Cv);
	assem4.push_data(ADV_param);
	assem4.push_vec(u);
	assem4.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem5(" Diff=data$2(#2);"
		  "V()+= comp (Base(#2))(i).Diff(i);"); 
	assem5.push_mi(mimv);
	assem5.push_mf(mf_Cv);
	assem5.push_mf(mf_coefv);
	assem5.push_data(Cv);
	assem5.push_data(DIFF_param);
	assem5.push_vec(a);
	assem5.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem6("l1=data$2(#2); "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l1(j);"); 
	assem6.push_mi(mimv);
	assem6.push_mf(mf_Cv);
	assem6.push_mf(mf_coefvi[i]);
	assem6.push_mf(mf_coefv);
	assem6.push_data(Cv);
	assem6.push_data(param.lambdax(i));
	assem6.push_data(param.lambday(i));
	assem6.push_data(param.lambdaz(i));
	assem6.push_data(DIFF_param);
	assem6.push_vec(l1);
	assem6.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem7("l2=data$3(#2); "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l2(j);"); 
	assem7.push_mi(mimv);
	assem7.push_mf(mf_Cv);
	assem7.push_mf(mf_coefvi[i]);
	assem7.push_mf(mf_coefv);
	assem7.push_data(Cv);
	assem7.push_data(param.lambdax(i));
	assem7.push_data(param.lambday(i));
	assem7.push_data(param.lambdaz(i));
	assem7.push_data(DIFF_param);
	assem7.push_vec(l1);
	assem7.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem9("l3=data$4(#2);  "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l3(j);"); 
	assem9.push_mi(mimv);
	assem9.push_mf(mf_Cv);
	assem9.push_mf(mf_coefvi[i]);
	assem9.push_mf(mf_coefv);
	assem9.push_data(Cv);
	assem9.push_data(param.lambdax(i));
	assem9.push_data(param.lambday(i));
	assem9.push_data(param.lambdaz(i));
	assem9.push_data(DIFF_param);
	assem9.push_vec(l1);
	assem9.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem8("c=data$1(#1); "
		  "t=comp(Grad(#1));"
		  "V$1(3)+=t(i,:).c(i);"); 
	assem8.push_mi(mimv);
	assem8.push_mf(mf_Cv);
	assem8.push_mf(mf_coefvi[i]);
	assem8.push_mf(mf_coefv);
	assem8.push_data(Cv);
	assem8.push_data(param.lambdax(i));
	assem8.push_data(param.lambday(i));
	assem8.push_data(param.lambdaz(i));
	assem8.push_data(DIFF_param);
	assem8.push_vec(cv);
	assem8.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem1("c=data$1(#1); l1=data$2(#2); l2=data$3(#2); l3=data$4(#2);  Diff=data$5(#3);"
		  "t=comp(Base(#3).Grad(#1).Base(#2));"
		  "V$1()+=t(d,i,1,j).Diff(d).c(i).l1(j)+t(d,i,2,j).Diff(d).c(i).l2(j)+t(d,i,3,j).Diff(d).c(i).l3(j);"); 
	assem1.push_mi(mimv);
	assem1.push_mf(mf_Cv);
	assem1.push_mf(mf_coefvi[i]);
	assem1.push_mf(mf_coefv);
	assem1.push_data(Cv);
	assem1.push_data(param.lambdax(i));
	assem1.push_data(param.lambday(i));
	assem1.push_data(param.lambdaz(i));
	assem1.push_data(DIFF_param);
	assem1.push_vec(Diff);
	assem1.assembly(Jv_transp[j].rg);

	generic_assembly 
	assem2("c=data$1(#1); Adv=data$2(#2);"
		  "V()+= comp (Base(#1).Base(#2))(i,j).c(i).Adv(j);"); 
	assem2.push_mi(mimv);
	assem2.push_mf(mf_Cv);
	assem2.push_mf(mf_Uvi[i]);
	assem2.push_data(Cv);
	assem2.push_data(ADV_param);
	assem2.push_vec(Adv);
	assem2.assembly(Jv_transp[j].rg);

cout<<"MBD partial calcolata in mass_balance = " << Diff[0] <<endl;
cout<<"MBA partial calcolata in mass_balance = " << Adv[0] <<endl;
cout<<"c = "<< c[0]<< "; dx(c) = "<< cv[0]<< "; dy(c) = "<< cv[1]<< "; dz(c) = "<< cv[2] <<"; "<<endl;
cout<<"diffusion = "<< a[0]<< "; advection = "<< u[0]<< "; l1 "<< l1[0]<< "; l2 = "<< l2[0]<< "; l3 = "<< l3[0] <<"; "<<endl;

vector_type Uv_basis( mf_Uvi[i].nb_dof()); gmm::clear(Uv_basis);
vector_type Cv_basis(dof_transp.Cv()); gmm::clear(Cv_basis);

generic_assembly
assem10( "V(#1)+=comp(Base(#1))(:);");
assem10.push_mi(mimv);
assem10.push_mf(mf_Cv);
assem10.push_vec(Cv_basis);
assem10.assembly(Jv_transp[j].rg);

generic_assembly
assem11("V(#1)+=comp(Base(#1))(:);");
assem11.push_mi(mimv);
assem11.push_mf(mf_Uvi[i]);
assem11.push_vec(Uv_basis);
assem11.assembly(Jv_transp[j].rg);


	string branch_suff = "";
	std::ostringstream convert;
	convert << i;
	branch_suff = convert.str();

	string junc_suff = "";
	std::ostringstream convert2;
	convert2 << j;
	junc_suff = convert2.str();

	std::ofstream outC("C_branch_"+branch_suff+"_junc_"+junc_suff+".txt");
		outC << gmm::col_vector(Cv_basis);
		outC.close();			

	std::ofstream outU("U_branch_"+branch_suff+"_junc_"+junc_suff+".txt");
		outU << gmm::col_vector(Uv_basis);
		outU.close();			


cout<<"MBD_partial = "<< DIFF[0]<< "; MBD_cum = "<< Jv_transp[j].MBD <<";   DIFF_param = "<<DIFF_param[0]<<endl;
cout<<"MBA_partial = "<< ADV[0]<< "; MBA_cum = "<< Jv_transp[j].MBA <<";   ADV_param = "<<ADV_param[0]<<endl;




}			
			



		} //end of junction loop
	} // end of branch loop	

*/


	// ADVECTION IN VESSELS 		
/*	size_type shift = 0;
	
	//Build tangent versors
	vector_type lambdax; // tangent versor: x component
	vector_type lambday; // tangent versor: y component
	vector_type lambdaz; // tangent versor: z component
	std::ifstream ifs(descr.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr.MESH_FILEV);
	asm_tangent_versor(ifs, lambdax, lambday, lambdaz);
	ifs.close();
	
	// Cycle over branches
	for (size_type i=0; i<nb_branches; ++i){
		
		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		//Velocities in the i-th branch
 		vector_type Uvi( mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
		gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);
		//Tangent versors of the i-th branch
		vector_type lambdax_K, lambday_K, lambdaz_K;
		for(size_type j=0; j<mf_coefvi[i].nb_dof(); j++)
		{
			lambdax_K.emplace_back(lambdax[i]);
			lambday_K.emplace_back(lambday[i]);
			lambdaz_K.emplace_back(lambdaz[i]);
		}
	
		//Build Bv
		asm_advection_network(Bv, mimv, mf_Cv, mf_coefvi[i], mf_Uvi[i], Uvi, lambdax_K, lambday_K, lambdaz_K, meshv.region(i) );
	}
*/


/*! Build the flux term at the junction in order to chek the mass balance
    @f$ M=\int_{\mathcal{E}_{MIX}} \beta~u~v~d\sigma@f$ and
    @f$ F=\int_{\mathcal{E}_{MIX}} \beta~c0~v~d\sigma@f$
 */
/*!
	@param F        BC contribution to rhs
	@param M        BC contribution to mass matrix
	@param mim      The integration method to be used
	@param mf_c     The finite element method for the concentration @f$u@f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of network boundary points
	@param beta     The beta value for mix condition @f$p_0@f$
	@ingroup asm
 */ 
template<typename VEC, typename VEC2>
void
asm_mass_fluxes
	(VEC & D, VEC & A,
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_lambda,
	 const mesh_fem & mf_dataD,
	 const mesh_fem & mf_dataA,
	 const VEC2 & C,
	 const VEC2 & lambdax, const VEC2 & lambday, const VEC2 & lambdaz,
	 const VEC2 & D_param,
	 const VEC2 & A_param,
	 const mesh_region & rg	
	) 
{ 
	GMM_ASSERT1(mf_c.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");

	//Build the exchange of mass at junction rg due to diffusion fluxes	

	generic_assembly 
	assem1("c=data$1(#1); l1=data$2(#2); l2=data$3(#2); l3=data$4(#2);  Diff=data$5(#3);"
		  "t=comp(Base(#3).Grad(#1).Base(#2));"
		  "V$1()+=t(d,i,1,j).Diff(d).c(i).l1(j)+t(d,i,2,j).Diff(d).c(i).l2(j)+t(d,i,3,j).Diff(d).c(i).l3(j);"); 
	assem1.push_mi(mim);
	assem1.push_mf(mf_c);
	assem1.push_mf(mf_lambda);
	assem1.push_mf(mf_dataD);
	assem1.push_data(C);
	assem1.push_data(lambdax);
	assem1.push_data(lambday);
	assem1.push_data(lambdaz);
	assem1.push_data(D_param);
	assem1.push_vec(D);
	assem1.assembly(rg);


	generic_assembly 
	assem2("c=data$1(#1); Adv=data$2(#2);"
		  "V()+= comp (Base(#1).Base(#2))(i,j).c(i).Adv(j);"); 
	assem2.push_mi(mim);
	assem2.push_mf(mf_c);
	assem2.push_mf(mf_dataA);
	assem2.push_data(C);
	assem2.push_data(A_param);
	assem2.push_vec(A);
	assem2.assembly(rg);

std::vector<scalar_type> c(1);
std::vector<scalar_type> cv(3);
std::vector<scalar_type> u(1);
std::vector<scalar_type> a(1);
std::vector<scalar_type> l1(1);
std::vector<scalar_type> l2(1);
std::vector<scalar_type> l3(1);


	generic_assembly 
	assem3("c=data$1(#1);"
		  "V()+= comp (Base(#1))(i).c(i);"); 
	assem3.push_mi(mim);
	assem3.push_mf(mf_c);
	assem3.push_mf(mf_dataA);
	assem3.push_data(C);
	assem3.push_data(A_param);
	assem3.push_vec(c);
	assem3.assembly(rg);

	generic_assembly 
	assem4(" Adv=data$2(#2);"
		  "V()+= comp (Base(#2))(i).Adv(i);"); 
	assem4.push_mi(mim);
	assem4.push_mf(mf_c);
	assem4.push_mf(mf_dataA);
	assem4.push_data(C);
	assem4.push_data(A_param);
	assem4.push_vec(u);
	assem4.assembly(rg);

	generic_assembly 
	assem5(" Diff=data$2(#2);"
		  "V()+= comp (Base(#2))(i).Diff(i);"); 
	assem5.push_mi(mim);
	assem5.push_mf(mf_c);
	assem5.push_mf(mf_dataD);
	assem5.push_data(C);
	assem5.push_data(D_param);
	assem5.push_vec(a);
	assem5.assembly(rg);

	generic_assembly 
	assem6("l1=data$2(#2); "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l1(j);"); 
	assem6.push_mi(mim);
	assem6.push_mf(mf_c);
	assem6.push_mf(mf_lambda);
	assem6.push_mf(mf_dataD);
	assem6.push_data(C);
	assem6.push_data(lambdax);
	assem6.push_data(lambday);
	assem6.push_data(lambdaz);
	assem6.push_data(D_param);
	assem6.push_vec(l1);
	assem6.assembly(rg);

	generic_assembly 
	assem7("l2=data$3(#2); "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l2(j);"); 
	assem7.push_mi(mim);
	assem7.push_mf(mf_c);
	assem7.push_mf(mf_lambda);
	assem7.push_mf(mf_dataD);
	assem7.push_data(C);
	assem7.push_data(lambdax);
	assem7.push_data(lambday);
	assem7.push_data(lambdaz);
	assem7.push_data(D_param);
	assem7.push_vec(l2);
	assem7.assembly(rg);

	generic_assembly 
	assem9("l3=data$4(#2);  "
		  "t=comp(Base(#2));"
		  "V$1()+=t(j).l3(j);"); 
	assem9.push_mi(mim);
	assem9.push_mf(mf_c);
	assem9.push_mf(mf_lambda);
	assem9.push_mf(mf_dataD);
	assem9.push_data(C);
	assem9.push_data(lambdax);
	assem9.push_data(lambday);
	assem9.push_data(lambdaz);
	assem9.push_data(D_param);
	assem9.push_vec(l3);
	assem9.assembly(rg);

	generic_assembly 
	assem8("c=data$1(#1); "
		  "t=comp(Grad(#1));"
		  "V$1(3)+=t(i,:).c(i);"); 
	assem8.push_mi(mim);
	assem8.push_mf(mf_c);
	assem8.push_mf(mf_lambda);
	assem8.push_mf(mf_dataD);
	assem8.push_data(C);
	assem8.push_data(lambdax);
	assem8.push_data(lambday);
	assem8.push_data(lambdaz);
	assem8.push_data(D_param);
	assem8.push_vec(cv);
	assem8.assembly(rg);


cout<<"c = "<< c[0]<< "; dx(c) = "<< cv[0]<< "; dy(c) = "<< cv[1]<< "; dz(c) = "<< cv[2] <<"; "<<endl;
cout<<"diffusion = "<< a[0]<< "; advection = "<< u[0]<< "; l1 "<< l1[0]<< "; l2 = "<< l2[0]<< "; l3 = "<< l3[0] <<"; "<<endl;
cout << "                                ---------- "   << endl;



}



//! Compute mass matrix with three parameters, that is:
//! @f$ M = \int_{\Omega} p_1~/left(~p_2~+~p_3/rigt)~u~v~dx @f$ and
/*!
	@param M		Computed mass matrix
	@param mim		The integration metod used
	@param mf_c		The finite element method for @f$ u @f$ and @f$ v @f$
	@param mf_data1		The finite element method for parameter @f$ p_1 @f$
	@param mf_data2		The finite element method for parameter @f$ p_2 @f$
	@param mf_data3		The finite element method for parameter @f$ p_3 @f$
	@param PARAM1		First parameter @f$ p_1 @f$
	@param PARAM2		Second parameter @f$ p_2 @f$
	@param PARAM3		Third parameter @f$ p_3 @f$
	@param rg		The region where to integrate

*/
template<typename MAT, typename VEC1, typename VEC2, typename VEC3>
void 
asm_mass_matrix_three_param
	(MAT & M, 
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data1,
	 const mesh_fem & mf_data2,
	 const mesh_fem & mf_data3,
	 const VEC1 & PARAM1,
	 const VEC2 & PARAM2,
	 const VEC3 & PARAM3,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{

	generic_assembly 
	assem("param1=data$1(#1); param2=data$2(#2); param3=data$3(#3);"
		"M$1(#4,#4)+=comp(Base(#4).Base(#4).Base(#1).Base(#2))(:,:,i,j).param1(i).param2(j);" 			"M$1(#4,#4)+=comp(Base(#4).Base(#4).Base(#1).Base(#3))(:,:,i,j).param1(i).param3(j);");
	assem.push_mi(mim);
	assem.push_mf(mf_data1);
	assem.push_mf(mf_data2);
	assem.push_mf(mf_data3);
	assem.push_mf(mf_c);
	assem.push_data(PARAM1);
	assem.push_data(PARAM2);	
	assem.push_data(PARAM3);
	assem.push_mat(M);
	assem.assembly(rg);
}



/*!
	Build the exchange matrices
	@f$B_{tt} = 2\pi ~R~ \left( - A + B \right)~ \Pi^T_{tv} M_{vv} \bar{\Pi}_{tv}@f$,
	@f$B_{tv} = 2\pi ~R~ \left( - A - B \right)~ \Pi^T_{tv} M_{vv}@f$,
	@f$B_{vt} = 2\pi ~R~ \left( + A - B \right)~ M_{vv} \bar{\Pi}_{tv}@f$,
	@f$B_{vv} = 2\pi ~R~ \left( + A + B \right)~ M_{vv}@f$,
	where @f$ A @f$ is the oncotic term and @f$ B @f$ is the vessel permeability term.   
	If ALT_FORM==true we substitute @f${\Pi}_{tv}@f$ with @f$\bar{\Pi}_{tv}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_exchange_mat_transp
	(MAT & Btt, MAT & Btv, MAT & Bvt, MAT & Bvv, 
	 const getfem::mesh_im & mim,
	 const getfem::mesh_fem & mf_c, 
	 const getfem::mesh_fem & mf_coefv,
	 const getfem::mesh_fem & mf_pv,
	 const MAT & Mbar, const MAT & Mlin,
	 const VEC & R,
	 const VEC & ONCOTIC,
	 const VEC & PERM,
	 const bool ALT_FORM
	 ) 
{


	MAT Bvv_temp(mf_c.nb_dof(),mf_c.nb_dof()); gmm::clear(Bvv_temp);

	#ifdef M3D1D_VERBOSE_
	cout << "    Assembling Bvv ..." << endl;  
	#endif
	getfem::asm_mass_matrix_three_param(Bvv, mim, 
					    mf_c, mf_coefv, mf_pv, mf_coefv,
					gmm::scaled(R, +2*pi), 
					gmm::scaled(ONCOTIC, +1),
					gmm::scaled(PERM, +1) ); 
	#ifdef M3D1D_VERBOSE_
	cout << "    Assembling Bvt ..." << endl;
	#endif
	getfem::asm_mass_matrix_three_param(Bvv_temp, mim, 
					    mf_c, mf_coefv, mf_pv, mf_coefv,
					gmm::scaled(R, +2*pi), 
					gmm::scaled(ONCOTIC, +1),
					gmm::scaled(PERM, -1) ); 
	gmm::mult(Bvv_temp, Mbar, Bvt);
	gmm::clear(Bvv_temp);

	if (ALT_FORM){
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btv (alternative form) ..." << endl;
		#endif
		
		getfem::asm_mass_matrix_three_param(Bvv_temp, mim, 
					    mf_c, mf_coefv, mf_pv, mf_coefv,
					gmm::scaled(R, +2*pi), 
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, -1) ); 
		gmm::mult(gmm::transposed(Mbar), Bvv_temp, Btv); 
		gmm::clear(Bvv_temp);
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btt (alternative form) ..." << endl;
		#endif
		getfem::asm_mass_matrix_three_param(Bvv_temp, mim, 
					    mf_c, mf_coefv, mf_pv, mf_coefv,
					gmm::scaled(R, +2*pi), 
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, +1) ); 

		gmm::mult3(gmm::transposed(Mbar),Bvv_temp,Mbar, Btt);
		gmm::clear(Bvv_temp);
	}
	else{
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btv (alternative form) ..." << endl;
		#endif
		
		getfem::asm_mass_matrix_three_param(Bvv_temp, mim, 
					    mf_c, mf_coefv, mf_pv, mf_coefv,
					gmm::scaled(R, +2*pi), 
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, -1) ); 
		gmm::mult(gmm::transposed(Mlin), Bvv_temp, Btv); 
		gmm::clear(Bvv_temp);
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btt (alternative form) ..." << endl;
		#endif
		getfem::asm_mass_matrix_three_param(Bvv_temp, mim, 
					    mf_c, mf_coefv, mf_pv, mf_coefv,
					gmm::scaled(R, +2*pi), 
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, +1) ); 

		gmm::mult3(gmm::transposed(Mlin),Bvv_temp,Mbar, Btt);
		gmm::clear(Bvv_temp);
	}
	
} /* end of build_exchange_matrices */


/// 3D test rhs
//! Exact concentration for 3d test
double c_exact_3d_1(const bgeot::base_node & x){
	return (1-x[0])*(exp(1)-exp(x[1]));
}
double c_exact_3d_2(const bgeot::base_node & x){
	return (1-x[0])*(exp(1)-exp(x[1]))/exp(1);
}
double c_exact_3d_3(const bgeot::base_node & x){
	return (1-x[0])*(1-x[1]);
}
double c_exact_3d_4(const bgeot::base_node & x){
	return (1-x[0])*(1-x[1])*(1-x[2]);
}
double c_exact_3d_5(const bgeot::base_node & x){
	return (sin(2*pi*x[0]))*(sin(2*pi*x[1]))*(sin(2*pi*x[2]));
}
double c_exact_3d_6(const bgeot::base_node & x){
	return (1-x[0])*(1-x[1])*(1-x[2])*(1-x[0])*(1-x[1])*(1-x[2]);
}
	if (PARAM.int_value("TEST_RHS_TRANSP")!=0) {
		#ifdef M3D1D_VERBOSE_
		cout << "  ... with an exact concentration ... " << endl;
		#endif
	vector_type c_ex(dof_transp.Ct());
		if (PARAM.int_value("TEST_RHS_TRANSP")==1){
	cout <<"test rhs transp number: 1"<<endl;
	interpolation_function(mf_Ct, c_ex, c_exact_3d_1 );}
		if (PARAM.int_value("TEST_RHS_TRANSP")==2){
	cout <<"test rhs transp number: 2"<<endl;
	interpolation_function(mf_Ct, c_ex, c_exact_3d_2 );}
		if (PARAM.int_value("TEST_RHS_TRANSP")==3){
	cout <<"test rhs transp number: 3"<<endl;	
	interpolation_function(mf_Ct, c_ex, c_exact_3d_3 );}
		if (PARAM.int_value("TEST_RHS_TRANSP")==4){
	cout <<"test rhs transp number: 4"<<endl;
	interpolation_function(mf_Ct, c_ex, c_exact_3d_4 );}
		if (PARAM.int_value("TEST_RHS_TRANSP")==5){
	cout <<"test rhs transp number: 5"<<endl;
	interpolation_function(mf_Ct, c_ex, c_exact_3d_5 );}
	if (PARAM.int_value("TEST_RHS_TRANSP")==6){
	cout <<"test rhs transp number: 6"<<endl;
	interpolation_function(mf_Ct, c_ex, c_exact_3d_6 );}
	

	vector_type Ft(dof_transp.Ct());
	sparse_matrix_type Att(dof_transp.Ct(), dof_transp.Ct());

	gmm::add(	gmm::sub_matrix(AM_temp,
			gmm::sub_interval(0,dof_transp.Ct()),
			gmm::sub_interval(0,dof_transp.Ct()))
			, Att);
	gmm::scale(	gmm::sub_matrix(AM_temp,
			gmm::sub_interval(0,dof_transp.Ct()),
			gmm::sub_interval(0,dof_transp.Ct()))
			, 0.0);	
			
	gmm::add(	 gmm::sub_vector(FM_temp, gmm::sub_interval(0,dof_transp.Ct()))
			,Ft);	 
	gmm::scale(	 gmm::sub_vector(FM_temp, gmm::sub_interval(0,dof_transp.Ct()))
			,0.0);
	
		for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
				getfem::assembling_Dirichlet_condition(Att, Ft, mf_Ct, BCt_transp[bc].rg, c_ex);}

	gmm::add(Att, 
			gmm::sub_matrix(AM_temp,
					gmm::sub_interval(0,dof_transp.Ct()),
					gmm::sub_interval(0,dof_transp.Ct())));
	gmm::add(Ft, 
			gmm::sub_vector(FM_temp,
					gmm::sub_interval(0,dof_transp.Ct())));
	// De-allocate memory
	gmm::clear(Att);
	gmm::clear(Ft);

	}else{




/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                         A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   assembling3d1d_transp.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 -May 2018.
  @brief  Miscelleanous assembly routines for the 3D/1D coupling for transport problem.
 */


#ifndef M3D1D_ASSEMBLING_3D1D_TRANSP_HPP_
#define M3D1D_ASSEMBLING_3D1D_TRANSP_HPP_

#include <defines.hpp>
#include <utilities.hpp>
#include <utilities_transp.hpp>

namespace getfem {


//! Compute mass matrix with two parameters, that is:
//! @f$ M = \int_{\Omega} ~p_3~(~p_1~+~p_2 )~u~v~dx @f$ and
/*!
	@param M		Computed mass matrix
	@param mim		The integration metod used
	@param mf_c		The finite element method for @f$ u @f$ and @f$ v @f$
	@param mf_data1		The finite element method for parameter @f$ p_1 @f$
	@param mf_data2		The finite element method for parameter @f$ p_2 @f$
	@param mf_data3		The finite element method for parameter @f$ p_3 @f$
	@param PARAM1		First parameter @f$ p_1 @f$
	@param PARAM2		Second parameter @f$ p_2 @f$
	@param PARAM3		Third parameter @f$ p_3 @f$	
	@param rg		The region where to integrate

*/
template<typename MAT, typename VEC1, typename VEC2, typename VEC3>
void 
asm_mass_matrix_three_param
	(MAT & M, 
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data1,
	 const mesh_fem & mf_data2,
	 const mesh_fem & mf_data3,
	 const VEC1 & PARAM1,
	 const VEC2 & PARAM2,
	 const VEC3 & PARAM3,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{

	generic_assembly 
	assem("param1=data$1(#1); param2=data$2(#2); param3=data$3(#3);"
		"M$1(#4,#4)+=comp(Base(#4).Base(#4).Base(#1).Base(#3))(:,:,i,j).param1(i).param3(j);" 
		"M$1(#4,#4)+=comp(Base(#4).Base(#4).Base(#2).Base(#3))(:,:,i,j).param2(i).param3(j);");
	assem.push_mi(mim);
	assem.push_mf(mf_data1);
	assem.push_mf(mf_data2);
	assem.push_mf(mf_data3);
	assem.push_mf(mf_c);
	assem.push_data(PARAM1);
	assem.push_data(PARAM2);
	assem.push_data(PARAM3);
	assem.push_mat(M);
	assem.assembly(rg);
}



/*!
	Build the exchange matrices
	@f$B_{tt} = \left( - A + B \right)~ \Pi^T_{tv} M_{vv} \bar{\Pi}_{tv}@f$,    
	@f$B_{tv} = \left( - A - B \right)~ \Pi^T_{tv} M_{vv}@f$,    
	@f$B_{vt} = \left( + A - B \right)~ M_{vv} \bar{\Pi}_{tv}@f$,    
	@f$B_{vv} = \left( + A + B \right)~ M_{vv}@f$,    
	where @f$ A @f$ is the oncotic term and @f$ B @f$ is the vessel permeability term.   
	If ALT_FORM==true we substitute @f${\Pi}_{tv}@f$ with @f$\bar{\Pi}_{tv}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_exchange_mat_transp
	(MAT & Btt, MAT & Btv, MAT & Bvt, MAT & Bvv, 
	 const getfem::mesh_im & mim,
	 const getfem::mesh_fem & mf_c, 
	 const getfem::mesh_fem & mf_ONCOTIC,
	 const getfem::mesh_fem & mf_PERM,
	 const getfem::mesh_fem & mf_R,
	 const MAT & Mbar, const MAT & Mlin,
	 const VEC & ONCOTIC,
	 const VEC & PERM,
	 const VEC & R,
	 const bool ALT_FORM
	 ) 
{


	MAT Bvv_temp(mf_c.nb_dof(),mf_c.nb_dof()); gmm::clear(Bvv_temp);

	#ifdef M3D1D_VERBOSE_
	cout << "    Assembling Bvv ..." << endl;  
	#endif
	getfem::asm_mass_matrix_three_param(Bvv, mim, 
					    mf_c, mf_ONCOTIC, mf_PERM, mf_R,
					gmm::scaled(ONCOTIC, +1),
					gmm::scaled(PERM, +1),
					gmm::scaled(R, 2*pi) ); 
	#ifdef M3D1D_VERBOSE_
	cout << "    Assembling Bvt ..." << endl;
	#endif
	getfem::asm_mass_matrix_three_param(Bvv, mim, 
					    mf_c, mf_ONCOTIC, mf_PERM, mf_R,
					gmm::scaled(ONCOTIC, +1),
					gmm::scaled(PERM, -1),
					gmm::scaled(R, 2*pi) ); 
	gmm::mult(Bvv_temp, Mbar, Bvt);
	gmm::clear(Bvv_temp);

	if (ALT_FORM){
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btv (alternative form) ..." << endl;
		#endif
		
	getfem::asm_mass_matrix_three_param(Bvv, mim, 
					    mf_c, mf_ONCOTIC, mf_PERM, mf_R,
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, -1),
					gmm::scaled(R, 2*pi) ); 
		gmm::mult(gmm::transposed(Mbar), Bvv_temp, Btv); 
		gmm::clear(Bvv_temp);
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btt (alternative form) ..." << endl;
		#endif
	getfem::asm_mass_matrix_three_param(Bvv, mim, 
					    mf_c, mf_ONCOTIC, mf_PERM, mf_R,
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, +1),
					gmm::scaled(R, 2*pi) ); 

		gmm::mult3(gmm::transposed(Mbar),Bvv_temp,Mbar, Btt);
		gmm::clear(Bvv_temp);
	}
	else{
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btv (alternative form) ..." << endl;
		#endif
		
	getfem::asm_mass_matrix_three_param(Bvv, mim, 
					    mf_c, mf_ONCOTIC, mf_PERM, mf_R,
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, -1),
					gmm::scaled(R, 2*pi) ); 
		gmm::mult(gmm::transposed(Mlin), Bvv_temp, Btv); 
		gmm::clear(Bvv_temp);
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btt (alternative form) ..." << endl;
		#endif
	getfem::asm_mass_matrix_three_param(Bvv, mim, 
					    mf_c, mf_ONCOTIC, mf_PERM, mf_R,
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, -1),
					gmm::scaled(R, 2*pi) ); 

		gmm::mult3(gmm::transposed(Mlin),Bvv_temp,Mbar, Btt);
		gmm::clear(Bvv_temp);
	}
	
} /* end of build_exchange_matrices */



  //! Build a single tissue Dirichlet condition on the network (modify @f$ B_{vt} @f$)
  /* @f{equation}{
	B = 
	@f{bmatrix} {B_{tt}& B_{tv}\\
		     B_{vt}& B_{vv}  @f}
	}
  */
  //! @f$ BU = F @f$

  /*!

	@param B	The monolitic matrix of the system
	@param F	The right hand side vector of the system  
	@param mf1 	The finite element method on tissue
	@param mf2 	The finite element method on network
	@param boundary	The index of the boundary region to be buildt
	@param DIR	The vector containing the Dirichlet condition (should be of the same dimension of F)

     @ingroup asm
  */
 template<typename MATRM, typename VECT1, typename VECT2>
  void assembling_Dirichlet_condition_coupled_tissue
  (MATRM &B, VECT1 &F, const mesh_fem &mf1, const mesh_fem &mf2, size_type boundary,
   const VECT2 &DIR) {
   

    size_type Q1=mf1.get_qdim();
    size_type Q2=mf2.get_qdim();

    size_type nb_dof1=mf1.nb_dof();
    size_type nb_dof2=mf2.nb_dof();    
    
    GMM_ASSERT1(!(mf1.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
    GMM_ASSERT1(!(mf2.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
						
    dal::bit_vector nndof = mf1.basic_dof_on_region(boundary);
    pfem pf1;
    
    for (dal::bv_visitor cv(mf1.convex_index()); !cv.finished(); ++cv) {	 	//per tutti i convessi cv della mesh 1
    
      pf1 = mf1.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf1->dim());
      size_type nbd = pf1->nb_dof(cv);	      
      for (size_type i = 0; i < nbd; i++) {					//per tutti i dof i del convesso cv
	size_type dof1 = mf1.ind_basic_dof_of_element(cv)[i*Q1];				//trova l'indice delle colonne riferite all
	if (nndof.is_in(dof1) && pf1->dof_types()[i] == ldof) {			//se il dof i del convesso cv è in "boundary"
  
	  for (size_type j = nb_dof1; j < nb_dof1+ nb_dof2; j++) {				//allora per tutti i dof j della mesh 2
		for (size_type l = 0; l < Q1; ++l) {
			F[j] -= B(j, dof1+l) * DIR[dof1+l];
	    		B(j, dof1+l) =  0;
	    		B(dof1+l, j) =  0;
	    	}
	    }
	  } 
	}
     }
   } /* end of assembling_Dirichlet_condition_coupled_tissue*/

  //! Build a single vessel Dirichlet condition on the tissue (modify @f$ B_{tv} @f$)
  /* @f{equation}{
	B = 
	@f{bmatrix} {B_{tt} & B_{tv}\\ B_{vt} & B_{vv}  @f}
	}
  */
  //! @f$ BU = F @f$

  /*!

	@param B	The monolitic matrix of the system
	@param F	The right hand side vector of the system  
	@param mf1 	The finite element method on tissue
	@param mf2 	The finite element method on network
	@param boundary	The index of the boundary region to be buildt
	@param DIR	The vector containing the Dirichlet condition (should be of the same dimension of F)

     @ingroup asm
  */
 template<typename MATRM, typename VECT1, typename VECT2>
  void assembling_Dirichlet_condition_coupled_vessel
  (MATRM &B, VECT1 &F, const mesh_fem &mf1, const mesh_fem &mf2, size_type boundary,
   const VECT2 &DIR) {
   

    size_type Q1=mf1.get_qdim();
    size_type Q2=mf2.get_qdim();

    size_type nb_dof1=mf1.nb_dof();
    size_type nb_dof2=mf2.nb_dof();    
    
    GMM_ASSERT1(!(mf1.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
    GMM_ASSERT1(!(mf2.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
					
    dal::bit_vector nndof = mf2.basic_dof_on_region(boundary);
    pfem pf2;
    
    for (dal::bv_visitor cv(mf2.convex_index()); !cv.finished(); ++cv) {	 	//per tutti i convessi cv della mesh 1
    
      pf2 = mf2.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf2->dim());
      size_type nbd = pf2->nb_dof(cv);	      
      for (size_type i = 0; i < nbd; i++) {					//per tutti i dof i del convesso cv
	size_type dof2 = mf2.ind_basic_dof_of_element(cv)[i*Q2];				//trova l'indice delle colonne riferite all
	if (nndof.is_in(dof2) && pf2->dof_types()[i] == ldof) {			//se il dof i del convesso cv è in "boundary"
	  
	  for (size_type j = 0; j < nb_dof1; j++) {				//allora per tutti i dof j della mesh 2
		for (size_type l = 0; l < Q2; ++l) {
			F[j] -= B(j, nb_dof1 + dof2+l) * DIR[nb_dof1 + dof2+l];
	    		B(j, nb_dof1 + dof2+l) =  0;
	    		B(nb_dof1 + dof2+l, j) =  0;
	    	}
	    }
	  } 
	}
     }
   } /* end of assembling_Dirichlet_condition_coupled_vessel*/
   

  //! Build all the Dirichlet conditions on the coupling matrixes (modify both @f$ B_{tv} @f$ and @f$ B_{vt} @f$))
  /* @f{equation}{
	B = 
	@f{bmatrix} {B_{tt} & B_{tv}\\ B_{vt} & B_{vv}  @f}
	}
  */
  //! @f$ BU = F @f$

  /*!
	@param M		The monolitic matrix of the system
	@param F		The right hand side vector of the system  
	@param mf_ct 		The finite element method on tissue
	@param mf_cv 		The finite element method on network
	@param BC_tissue	List of nodes of the boundary regions in tissue
	@param BC_vessel	List of nodes of the boundary regions in network

     @ingroup asm
  */
template<typename MAT, typename VEC>
void
asm_coupled_bc_transp
	(MAT & M,
 	 VEC & F,
	 const mesh_fem & mf_ct,
	 const mesh_fem & mf_cv,
	 const std::vector<getfem::node> & BC_tissue,
	 const std::vector<getfem::node> & BC_vessel
	 )
{

	
	
	GMM_ASSERT1(mf_ct.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
	GMM_ASSERT1(mf_cv.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");


	//cycle over the tissue boundary nodes
	for (size_type bc=0; bc < BC_tissue.size(); ++bc) {
		GMM_ASSERT1(mf_ct.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC_tissue[bc].label=="DIR") { // Dirichlet BC
			VEC BC_temp(mf_ct.nb_dof(), BC_tissue[bc].value);
			getfem::assembling_Dirichlet_condition_coupled_tissue(M, F, mf_ct, mf_cv, BC_tissue[bc].rg, BC_temp);
			gmm::clear(BC_temp);				
		} 
	}
	
	//cycle over the vessels boundary nodes
	for (size_type bc=0; bc < BC_vessel.size(); ++bc) {
		GMM_ASSERT1(mf_cv.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC_vessel[bc].label=="DIR") { // Dirichlet BC
			VEC BC_temp(mf_cv.nb_dof(), BC_vessel[bc].value);
			getfem::assembling_Dirichlet_condition_coupled_vessel(M, F, mf_ct, mf_cv, BC_vessel[bc].rg, BC_temp);
			gmm::clear(BC_temp);				
		} 
	}


} /* end of asm_coupled_bc_transp */



} /* end of namespace */

#endif



	

	// bluid oncotic term
	vector_type ONCOTIC (dof.Pv());
	gmm::copy(gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
		  ONCOTIC);
	gmm::mult_add(gmm::scaled(Mbar,-1.0), 
		  gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut(), dof.Pt())),
		  ONCOTIC);


	scalar_type picoef=param_transp.sigma()*(param_transp.pi_v()-param_transp.pi_t());
        vector_type DeltaPi(dof.Pv(),picoef);
        gmm::add(gmm::scaled(DeltaPi,-1.0), ONCOTIC);	

	gmm::scale(ONCOTIC,0.5*(1.0-param_transp.sigma())*param.Q(0));
	
	// build permeability term
	vector_type PERM (dof.coefv(),param_transp.Y()[0]);

	//build exchange matrixes	
	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_Pv, mf_coefv, mf_coefv,  Mbar, Mlin, 
			ONCOTIC, PERM, param.R(), NEWFORM);
cout <<"R:  " <<param.R()<<endl;
cout <<"ONCOTIC:  " <<ONCOTIC<<endl;
cout <<"PERM:  " <<PERM<<endl;


