#ifndef __FV_H__
#define __FV_H__

//------------------------------------------------------------------------------
// Main problem class
//------------------------------------------------------------------------------

#include <vector>
#include <fstream>
#include <string>
#include "ext_constants.h"
#include "parameter.h"
#include "vec.h"
#include "primvar.h"
#include "flux.h"
#include "convar.h"
#include "entvar.h"
#include "misc_func.h"
#include "neural_network.h"

class FV
{
   public:
      FV(char* file)
      {
         param.input_file = file;
         param.read();
         SOLN_DIR = param.OUTPUT_PATH + "/SOLN_FILES_" + param.reconstruction;
         
         // Making solution directories
		 string commandline0 = "rm -rf " + SOLN_DIR;
		 system(commandline0.c_str());
		 string commandline1 = "mkdir -p " + SOLN_DIR;
		 system(commandline1.c_str());
		 
         glob_filename = SOLN_DIR + "/globals.dat";
		 glob_file.open (glob_filename.c_str());
		 
		 L1_err_filename = SOLN_DIR + "/l1_err.dat";
		 L2_err_filename = SOLN_DIR + "/l2_err.dat";
		 Linf_err_filename = SOLN_DIR + "/linf_err.dat";
		 
		 if(param.exact_available)
		 {
		    L1_err_file.open (L1_err_filename.c_str());
		    L2_err_file.open (L2_err_filename.c_str());
		    Linf_err_file.open (Linf_err_filename.c_str());
		 }
         
         if(param.OOC)
         {
           print_globals = false;
           print_soln = false;
           disp_log = false;
           L1_err.resize(param.mesh_levels);
           L2_err.resize(param.mesh_levels);
           Linf_err.resize(param.mesh_levels);
         }  
         else
         {
           print_globals = true; 
           print_soln = true; 
           disp_log = true;
         }  
           
      };
      void run ();
      ~FV()
      {
         glob_file.close();
         if(param.OOC)
		 {
		    L1_err_file.close();
		    L2_err_file.close();
		    Linf_err_file.close();
		 }
      };
   
   private:
      Neural_Network net;
      Parameter param;
      std::ofstream glob_file, L1_err_file,L2_err_file,Linf_err_file;
      std::string glob_filename,L1_err_filename,L2_err_filename,Linf_err_filename;
      void initialize ();
      void compute_dt ();
      void con_to_prim ();
      void prim_to_con ();
      void prim_to_ent ();
      void reconstruct (const ReconstructionScheme reconstruct_scheme,
                        const double h,
                        const std::vector< std::vector<double> > & state,
                        double (& rec_statel)[NVAR],
						double (& rec_stater) [NVAR]) const;
      void reconstruct_first (const std::vector< std::vector<double> > & state,
                              double (& rec_statel)[NVAR],
                              double (& rec_stater) [NVAR]) const ;
      
      void reconstruct_second_muscl (const std::vector< std::vector<double> > & state,
									 double (& rec_statel)[NVAR],
									 double (& rec_stater) [NVAR]) const;
      
      void reconstruct_muscl_limited(const ReconstructionScheme reconstruct_scheme,
       								 const double& h,
                                     const std::vector< std::vector<double> > & state,
									 double (& rec_statel)[NVAR],
									 double (& rec_stater) [NVAR]) const;  
	  void reconstruct_eno_pw(const int eno_level,
							const double& h,
							const std::vector< std::vector<double> > & state,
							double (& rec_statel)[NVAR],
							double (& rec_stater) [NVAR]) const;
	  void reconstruct_sp_weno(const double& h,
							const std::vector< std::vector<double> > & state,
							double (& rec_statel)[NVAR],
							double (& rec_stater) [NVAR]) const;
     void reconstruct_sp_weno_dl(const std::vector< std::vector<double> > & state,
							double (& rec_statel)[NVAR],
							double (& rec_stater) [NVAR]) const;															                     
      double limited_slope (const double& ul, 
                            const double& ur,
                            const ReconstructionScheme reconstruct_scheme,
                            const double& h) const;
      double Entropy (const PrimVar& state) const;
      double Enthalpy (const PrimVar& state) const;
      
      void kepec_flux(const std::vector<PrimVar> &P,
					const Vector &unit_normal, 
					Flux &flux) const;
	  void kep_flux(const std::vector<PrimVar> &P,
					const Vector &unit_normal, 
					Flux &flux) const;
      void kep4_flux(const std::vector<PrimVar> &P,
					const Vector &unit_normal, 
					Flux &flux) const; 
      void keps_flux(const std::vector<PrimVar> &P,
					const Vector &unit_normal, 
					Flux &flux) const;
      void keps4_flux(const std::vector<PrimVar> &P,
					const Vector &unit_normal, 
					Flux &flux) const;                            	
	  void roe_flux(const std::vector<PrimVar> &P,
					   const Vector &unit_normal,
					   Flux &flux) const;
      void roe_ec_flux(const std::vector<PrimVar> &P,
					   const Vector &unit_normal,
					   Flux &flux) const;						   					
	  void kepec4_flux(const std::vector<PrimVar> &P,
					const Vector &unit_normal, 
					Flux &flux) const;	
	  void roe_ec4_flux(const std::vector<PrimVar> &P,
					const Vector &unit_normal, 
					Flux &flux) const;		
	  void kep_tecno_roe_flux(const std::vector<PrimVar> &P,
								const std::vector<EntVar> &V,
								const Vector &unit_normal,
								Flux &flux) const;										
	  void kepes_tecno_roe_flux(const std::vector<PrimVar> &P,
								const std::vector<EntVar> &V,
								const Vector &unit_normal,
								Flux &flux) const;	
	  void kepes_tecno4_roe_flux(const std::vector<PrimVar> &P,
								const std::vector<EntVar> &V,
								const Vector &unit_normal,
								Flux &flux) const;			
	  void roe_tecno_roe_flux(const std::vector<PrimVar> &P,
					          const std::vector<EntVar> &V,
					          const Vector &unit_normal,
					          Flux &flux) const;															                         
      void kepes_tecno_rusanov_flux(const std::vector<PrimVar> &P,
									const std::vector<EntVar> &V,
									const Vector &unit_normal,
									Flux &flux) const;
	  void kepes_tecno4_rusanov_flux(const std::vector<PrimVar> &P,
									const std::vector<EntVar> &V,
									const Vector &unit_normal,
									Flux &flux) const;	
	  void roe_tecno4_roe_flux(const std::vector<PrimVar> &P,
					   const std::vector<EntVar> &V,
					   const Vector &unit_normal,Flux &flux) const;															
      void num_flux (const std::vector<PrimVar> &P,
					 const std::vector<EntVar> &V,
					 const Vector &unit_normal,
					 Flux &flux) const;              
      void compute_face_values ();
      void compute_residual (const int rk);
      void compute_inv_residual (const int rk);
      void compute_visc_residual (const int rk);
      void nodal_avg_and_vel_derivatives (int i, int j, int k,
                                        PrimVar &prim_avg, Vector &Txyz_vel,
                                        Vector &Tyxz_vel, Vector &Tzxy_vel) const;

      void nodal_vlux (int i, int j, int k,
                           Flux &xflux, Flux &yflux, Flux &zflux) const;
      void nodal_vlux_KE (int i, int j, int k,
                           Flux &xflux, Flux &yflux, Flux &zflux) const;   
      void nodal_vlux_ES (int i, int j, int k,
                           Flux &xflux, Flux &yflux, Flux &zflux) const;                                       
      void compute_residual_norm ();
      void update_solution (const int rk);                                        
      void output ();
      void apply_bc();
      void compute_global();
      void compute_error(const int &ref_level);
      void compute_OOC();
      void clear_data_structs();
      double Temperature(const PrimVar &prim) const;
      void find_node_weights();
      void find_vorticity();
      void load_network(const std::string net_file_name);
      //double enthalpy(const std::vector<double>& prim) const;
      //double viscosity (const double T) const;
    
      int iter;
      double time;
      double dt,dx,dy,dz;
      std::vector< std::vector< std::vector<Vector> > >cell_cc;
      std::vector< std::vector< std::vector<double> > >node_w;
      std::vector< std::vector< std::vector<PrimVar> > >primitive;
      std::vector< std::vector< std::vector<PrimVar> > >primitive_old;
      std::vector< std::vector< std::vector<Flux> > >residual;
      std::vector< std::vector< std::vector<Flux> > >residual2;  // FOR RK4
      std::vector< std::vector< std::vector<ConVar> > >conserved;
      std::vector< std::vector< std::vector<ConVar> > >conserved_old;
      std::vector< std::vector< std::vector<EntVar> > >entropy_var;
      std::vector< std::vector< std::vector<Vector> > >vorticity;
      Flux res_norm;
      Flux res_norm0;
      double total_ent0;
      
      bool print_globals,disp_log,print_soln;
      
      int Nx_cell, Ny_cell, Nz_cell;
      int rec_ng_cell, ng_cell;
      int counter;
      
      std::vector<PrimVar> L1_err,L2_err,Linf_err;
      
      std::string SOLN_DIR;
      
      double Cp;
      
      //Evaluated at cell centers
      double total_ent;
      double KE;
      double KE_diss_pr;
   
      //Evaluated at nodes
      double Enstrophy;
      double KE_diss_stress;
};

#endif
