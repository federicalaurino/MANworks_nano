/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                         A.Y. 2015-2016
                  
                Copyright (C) 2017 Stefano Brambilla
======================================================================*/
/*! 
  @file   main.cpp  
  @author Stefano Brambilla <s.brambilla93@gmail.com>   
  @date   September 2016 
  @brief  Main program for test simulations.    
  @details
    We solve the coupled 3D/1D problem of fluid exchange between a 1D 
    network \Lambda and the 3D interstitial tissue \Omega
    
    *****************************************
      Benchmark : Y-shaped bifurcation network 
      Mixed finite elements approximation
      Monolithic resolution by SuperLU 3.0
    *****************************************
    
	See Section "Code verification: test-cases"
 */   
 	
 	#define M3D1D_VERBOSE_
#include <iostream>
#include <problem3d1d.hpp>  
#include <transport3d1d.hpp> 
#include <problemHT.hpp> 

//! main program
int main(int argc, char *argv[])   
{

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
 
	try {   
		// Declare a new problem 
		getfem::transport3d1d p; 
		  
		/// fluid problem: velocity field and pressure
		
		// Initialize the problem
		p.init_fluid(argc, argv);
		// Build the monolithic system 
		p.assembly_fluid();
		// Solve the problem 
		if (!p.solve_fluid()) GMM_ASSERT1(false, "solve procedure has failed");
		// Save results in .vtk format
		p.export_vtk_fluid();
		// Display some global results: mean pressures, total flow rate
		
		std::cout << "--- FINAL RESULTS -------------------------" << std::endl; 
		std::cout << "  Pt average            = " << p.mean_pt()   << std::endl;
		std::cout << "  Pv average            = " << p.mean_pv()   << std::endl;
		std::cout << "  Network-to-Tissue TFR = " << p.flow_rate() << std::endl;
		std::cout << "-------------------------------------------" << std::endl; 	
		
		
		     
		//transport problem: concentration  
		
		//initialize 
		p.init_transp(argc, argv);
		//assemble        
		p.assembly_transp();    
		//solve     
		if (!p.solve_transp()) GMM_ASSERT1(false, "solve procedure has failed");  // the export is in the solve at each time step 
		//check the mass balance at each junction
		p.mass_balance();
				      
		   	
		}  
      
	GMM_STANDARD_CATCH_ERROR;     
		 
	   
    
		    
		   
	return 0;    
	   
} /* end of main program */   
    
   
  

