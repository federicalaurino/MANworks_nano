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
      Benchmark : single-vessel network 
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

		p.model_error(argc,argv);

		/*
		p.model_error_A0(argc,argv);
		p.model_error_A1(argc,argv);
		p.model_error_A2(argc,argv);
		p.model_error_A3(argc,argv);
		*/

/*
// A1: U=U(s)

	//Primal		  
		//initialize 
		p.init_U1(argc, argv);
		//assemble        
		p.assembly_U1();    
		//solve     
		if (!p.solve_U1()) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		p.export_U1();

	//Dual
		//initialize 
		p.init_Z1(argc, argv);
		//assemble        
		p.assembly_Z1();    
		//solve     
		if (!p.solve_Z1()) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		p.export_Z1();

	//Model error
		//compute error
		p.compute_error_1();
		   
	
//A2: Omega+ == Omega

	//Primal		  
		//initialize 
		p.init_U2(argc, argv);
		//assemble        
		p.assembly_U2();    
		//solve     
		if (!p.solve_U2()) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		p.export_U2();

	//Dual
		//initialize 
		p.init_Z2(argc, argv);
		//assemble        
		p.assembly_Z2();    
		//solve     
		if (!p.solve_Z2()) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		p.export_Z2();

	//Model error
		//compute error
		p.compute_error_2();
		

//A3: neglect fluctuations

	//Primal		  
		//initialize 
		p.init_U3(argc, argv);
		//assemble        
		p.assembly_U3();    
		//solve     
		if (!p.solve_U3()) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		p.export_U3();

	//Dual
		//initialize 
		p.init_Z3(argc, argv);
		//assemble        
		p.assembly_Z3();    
		//solve     
		if (!p.solve_Z3()) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		p.export_Z3();

	//Model error
		//compute error
		p.compute_error_3();
*/

}  
      
	GMM_STANDARD_CATCH_ERROR;     
		 
	   
    
		    
		   
	return 0;    
	   
} /* end of main program */   
    
   
  

