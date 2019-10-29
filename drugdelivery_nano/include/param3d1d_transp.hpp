/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   param3d1d_transp.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the aux class for physical parameters.
 */
 
#ifndef M3D1D_PARAM3D1D_TRANSP_HPP_
#define M3D1D_PARAM3D1D_TRANSP_HPP_

#include <mesh1d.hpp>    // import_network_radius
#include <utilities.hpp> // compute_radius

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
/*!
	\ingroup input
 */
struct param3d1d_transp {

	// Dimensional physical parameters (microcirc applications)
	//Diffusivity in the tissue [m^2/s]
	scalar_type Dt_;
	//Diffusivity in the vessels [m^2/s]
	scalar_type Dv_;
	//rate of metabolization [1/s]
	scalar_type m_;
	//Permeability of the vessel wall [m/s]
	scalar_type Perm_;
	//hydraulic conductivity of the lymphatic wall [s * m^2/kg]
	scalar_type Lp_LF_;
	// surface area of lymphatic vessels per unit volume of tissue [1/m]
	scalar_type SV_;
	// characteristic velocity
	scalar_type U_;
	// characteristic length
	scalar_type d_;
	
	
	// Dimensionless physical parameters (test-cases)
	//Inverse of Peclet number for tissue
	vector_type At_;
	//Inverse of Peclet number for vessel
	vector_type Av_;
	//Damkohler number (metabolism VS diffusion)
	vector_type Dalpha_;
	//Magnitude of leakage from the capillary bed
	vector_type Y_;
	//lymphatic drainage
	vector_type Q_pl_;
	//density of the fluid [kg/m^3] 
	scalar_type rho_;

	//Nanoparticle parameters
	//ligand surface density [1 / m^2]
	scalar_type m_l_;
	//receptor surface density [1 / m^2]
	scalar_type m_r_;
	//affinity constant ligand-receptor [m^2]
	scalar_type Ka_;
	//max distance for a ligand-receptor bond [m]
	scalar_type h_0_;
	// chord at h_0 [m]
	scalar_type r_0_;
	//beta_nano= 6 F \lambda /k_B T  [s^2 /Kg m]
	scalar_type beta_nano_;
	//diameter of the nanoparticle [m]
	scalar_type dp_;
	//Probability of adhesion [-]
	scalar_type Pa_;
	//Ligand percentage
	scalar_type rho_l_;

	//maximum value for the surface adhesive density  [1/m^2]
	scalar_type psi_max_;

	
	//time parameters 
	// simulation time length
	scalar_type T_;
	// time step
	scalar_type dt_;
	// initial concentration in tissue
	scalar_type C0t_;
	// initial concentration in network
	scalar_type C0v_;
	// injection time 
	scalar_type inj_time_;


	//Power law saturation parameter [s^b] (J_coef)
	scalar_type ups_;
	//Power law exponent [-] (J_coef)
	scalar_type b_;
	//total drug concentration in NP [1/m^3] (J_coef)
	scalar_type cNP_;
	//NP volume [m^3] (J_coef)
	scalar_type VNP_;
	

	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Finite Element Method for tissue data
	getfem::mesh_fem mf_datat_;
	//! Finite Element Method for vessel data
	getfem::mesh_fem mf_datav_;
	// Methods
	//! Build the arrays of dimensionless parameters
	void build(ftool::md_param & fname,
			const getfem::mesh_fem & mf_datat,
			 const getfem::mesh_fem & mf_datav
			) 
	{
		FILE_ = fname;
		mf_datat_ = mf_datat;
		mf_datav_ = mf_datav;
		size_type dof_datat = mf_datat_.nb_dof();
		size_type dof_datav = mf_datav_.nb_dof();
		 
//		bool IMPORT_RADIUS = FILE_.int_value("IMPORT_RADIUS");
		bool NONDIM_PARAM  = FILE_.int_value("TEST_PARAM");
		bool EXPORT_PARAM  = FILE_.int_value("EXPORT_PARAM");
		bool PA_INPUT  = FILE_.int_value("PA_INPUT");
		bool SATURATION  = FILE_.int_value("SATURATION");
		bool TIME_INJECTION  = FILE_.int_value("TIME_INJECTION");
				
		
		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling dimensionless parameters Dt, Dv, Dalpha, Q_pl ... "   << endl;
		#endif
		if (NONDIM_PARAM) {
			// Import dimensionless params from FILE_
			scalar_type Atval = FILE_.real_value("At"); 
			scalar_type Avval = FILE_.real_value("Av"); 
			scalar_type Dalphaval = FILE_.real_value("D_alpha"); 
			scalar_type Yval = FILE_.real_value("Y"); 
			scalar_type Q_plval  = FILE_.real_value("Q_pl"); 
			// Fill the data arrays
			 At_.assign(dof_datat,  Atval);
			 Av_.assign(dof_datav,  Avval);
			 Dalpha_.assign(dof_datat,  Dalphaval);
			 Y_.assign(dof_datav,  Yval);
			 Q_pl_.assign(dof_datat,  Q_plval);

			 //Import nanoparticles parameters
						//NANO
			m_l_   = 	FILE_.real_value("m_l","ligand surface density [1 / m^2]");
			m_r_   = 	FILE_.real_value("m_r","receptor surface density [1 / m^2]");
			Ka_    = 	FILE_.real_value("Ka","affinity constant ligand-receptor [m^2]");
			h_0_ =		FILE_.real_value("h_0","max distance for a ligand-receptor bond [m]");
			r_0_	=	FILE_.real_value("r_0","chord at h_0 [m]");
			beta_nano_=  	FILE_.real_value("beta_nano","beta_nano= 6 F lambda /k_B T  [s^2 /Kg m]");
			dp_   = FILE_.real_value("dp","diameter of the nanoparticle [m]");
			rho_  = 	FILE_.real_value("rho", "Density of the fluid (blood) [kg/m^3] ");
			if (PA_INPUT){
			rho_l_   = 	FILE_.real_value("rho_l","ligand percentage [-]");
			}
			if (SATURATION){
			psi_max_   = 	FILE_.real_value("psi_max","maximum value for the surface adhesive density  [1/m^2]");
			}
			if(TIME_INJECTION){
			inj_time_= FILE_.real_value("inj_time", "injection time [s]");
			}

			//J_coef
			ups_= FILE_.real_value("ups", "Power law saturation parameter [s^b]"); 
			b_= FILE_.real_value("b", "Power law exponent [-]");
 			cNP_= FILE_.real_value("cNP", "total drug concentration in NP [1/m^3]"); 
			VNP_= FILE_.real_value("VNP", "NP volume [m^3]"); 


			T_   = FILE_.real_value("T","Simulation time length [s]");
			dt_   = FILE_.real_value("dt","Time step [s]");			
			C0t_   = FILE_.real_value("C0t","Initial concentration in tissue []");	
			C0v_   = FILE_.real_value("C0v","Initial concentration in network []");			
		} 
		else { 
			// Import dimensional params from FILE_
			scalar_type P_  = FILE_.real_value("P", "average interstitial pressure [Pa]"); 
			U_  = FILE_.real_value("U", "characteristic flow speed in the capillary bed [m/s]"); 
			d_  = FILE_.real_value("d", "characteristic length of the problem [m]"); 
			scalar_type k_  = FILE_.real_value("k", "permeability of the interstitium [m^2]"); 
			scalar_type Lp_ = FILE_.real_value("Lp", "Hydraulic conductivity of the capillary walls [m^2 s/kg]"); 
			
			Dt_   = FILE_.real_value("Dt","Diffusivity in the tissue [m^2/s]");
			Dv_   = FILE_.real_value("Dv","Diffusivity in the vessels [m^2/s]");
			m_    = FILE_.real_value("m","rate of metabolization [1/s]");
			Perm_ = FILE_.real_value("Perm","Permeability of the capillary walls [m/s]");
			Lp_LF_ = FILE_.real_value("Lp_LF","hydraulic conductivity of the lymphatic wall [s * m^2/kg]");
			SV_ = FILE_.real_value("SV","surface area of lymphatic vessels per unit volume of tissue [1/m]");
			rho_  = 	FILE_.real_value("rho", "Density of the fluid (blood) [kg/m^3] ");

			//Import nanoparticles parameters
						//NANO
			m_l_   = 	FILE_.real_value("m_l","ligand surface density [1 / m^2]");
			m_r_   = 	FILE_.real_value("m_r","receptor surface density [1 / m^2]");
			Ka_    = 	FILE_.real_value("Ka","affinity constant ligand-receptor [m^2]");
			h_0_ =		FILE_.real_value("h_0","max distance for a ligand-receptor bond [m]");
			r_0_	=	FILE_.real_value("r_0","chord at h_0 [m]");
			beta_nano_=  	FILE_.real_value("beta_nano","beta_nano= 6 F lambda /k_B T  [s^2 /Kg m]");
			dp_   = FILE_.real_value("dp","diameter of the nanoparticle [m]");
			if (PA_INPUT){
			rho_l_   = 	FILE_.real_value("rho_l","ligand percentage [-]");
			}
			if (SATURATION){
			psi_max_   = 	FILE_.real_value("psi_max","maximum value for the surface adhesive density  [1/m^2]");
			}
			if(TIME_INJECTION){
			inj_time_= FILE_.real_value("inj_time", "injection time [s]");
			}

			//J_coef
			ups_= FILE_.real_value("ups", "Power law saturation parameter [s^b]"); 
			b_= FILE_.real_value("b", "Power law exponent [-]");
 			cNP_= FILE_.real_value("cNP", "total drug concentration in NP [1/m^3]"); 
			VNP_= FILE_.real_value("VNP", "NP volume [m^3]"); 


			T_   = FILE_.real_value("T","Simulation time length [s]");
			dt_   = FILE_.real_value("dt","Time step [s]");	
			C0t_   = FILE_.real_value("C0t","Initial concentration in tissue []");	
			C0v_   = FILE_.real_value("C0v","Initial concentration in network []");				
			// Compute the dimentionless params
			At_.assign(dof_datat, Dt_/d_/U_);
			Av_.assign(dof_datav, Dv_/d_/U_);
			Dalpha_.assign(dof_datat, m_/U_*d_);
			Y_.assign(dof_datav, Perm_/U_);
			Q_pl_.assign(dof_datat,Lp_LF_*SV_*P_*d_/U_);
						
		}
		// Check values
		GMM_ASSERT1(At_[0] != 0, "wrong tissue diffusivity (At>0 required)"); 
		GMM_ASSERT1(Av_[0] != 0, "wrong vessel bed diffusivity (Av>0 required)");
		//if (Q_[0] == 0) cout << "Warning: uncoupled problem (Q=0)" << endl;
		


	}
	//! Get the radius at a given dof
	//inline scalar_type R  (size_type i) { return R_[i];  } const
	//! Get the tissue diffusivity at a given dof
	inline scalar_type At (size_type i) { return At_[i]; } const
	//! Get the vessel diffusivity at a given dof
	inline scalar_type Av (size_type i) { return Av_[i]; } const
	//! Get the linphatic drainage at a given dof
	inline scalar_type Q_pl  (size_type i) { return Q_pl_[i];  } const
	//! Get the Dahmkholer number at a given dof
	inline scalar_type Dalpha  (size_type i) { return Dalpha_[i];  } const
	//! Get the leakage of the capillary bed at a given dof
	inline scalar_type Y  (size_type i) { return Y_[i];  } const
	//! Get the blood density
	inline scalar_type rho () { return rho_;  } const
	//! Get the simulation time length
	inline scalar_type T  () { return T_;  } const
	//! Get the time step
	inline scalar_type dt  () { return dt_;  } const
	//! Get the sinitial concentration in tissue
	inline scalar_type C0t  () { return C0t_;  } const
	//! Get the sinitial concentration in network
	inline scalar_type C0v  () { return C0v_;  } const
	//! Get the characteristic velocity
	inline scalar_type Uadim  () { return U_;  } const
	//! Get the characteristic lenght
	inline scalar_type dadim  () { return d_;  } const

	//Import Nanoparticle parameters
	//! Get the ligand surface density
	inline scalar_type & m_l () { return m_l_;  } const
	//! Get the receptor surface density
	inline scalar_type  & m_r  () { return m_r_;  } const
	//! Get the affinity constant ligand-receptor
	inline scalar_type  & Ka  () { return Ka_;  } const
	//! Get the max distance for a ligand-receptor bond
	inline scalar_type & h_0 () { return h_0_;  } const
	//! Get the chord at h_0
	inline scalar_type  & r_0  () { return r_0_;  } const
	//! Get the /beta_nano= 6 F \lambda /k_B T
	inline scalar_type  & beta_nano  () { return beta_nano_;  } const
	//! Get the diameter of nanoparticle
	inline scalar_type dp  () { return dp_;  } const
	//! Get the Probability of adhesion
	inline scalar_type Pa  () { return Pa_;  } const	
	//! Get the ligand percentage
	inline scalar_type rho_l  () { return rho_l_;  } const
	//! Get the maximum value for the surface adhesive density
	inline scalar_type psi_max  () { return psi_max_;  } const
	//! Get the time of injection
	inline scalar_type inj_time  () { return inj_time_;  } const

	//! Get Power law saturation parameter [s^b]		
	inline scalar_type  & ups  () { return ups_;  } const
	//! Get the Power law exponent [-]		
	inline scalar_type  & b  () { return b_;  } const
	//! Get the total drug concentration in NP [1/m^3]	
	inline scalar_type  & cNP  () { return cNP_;  } const
	//! Get the NP volume [m^3]
	inline scalar_type  & VNP  () { return VNP_;  } const
	
	//! Get the radius
	//vector_type & R (void) { return R_; }
	//! Get the vessel wall permeabilities
	vector_type & Q_pl (void) { return Q_pl_; }
	//! Get the Dahmkholer number 
	vector_type & Dalpha (void) { return Dalpha_; }
	//! Get the tissue diffusivity 
	vector_type & At (void) { return At_; }
	//! Get the vessel diffusivity 
	vector_type & Av (void) { return Av_; }
	//! Get the leakage of the capillary bed
	vector_type & Y  (void) { return Y_;  } const
	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const param3d1d_transp & param
		)
	{ 
		out << "--- PHYSICAL PARAMS ------" << endl;
		out << "  At     : "                << param.At_[0] << endl; 
		out << "  Av : "                << param.Av_[0] << endl; 
		out << "  Y      : "                << param.Y_[0] << endl; 
		out << "  Q_pl : "                << param.Q_pl_[0] << endl; 
		out << "  D_alpha : "                << param.Dalpha_[0] << endl; 
		out << "  T : "                << param.T_ << endl; 
		out << "  dt : "                << param.dt_ << endl; 
		out << "--------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
