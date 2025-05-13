#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <cassert>
#include <cstdlib>
#include "parameter.h"

using namespace std;

extern map<string,double> constants;

//------------------------------------------------------------------------------
// Read parameters from file
//------------------------------------------------------------------------------
void Parameter::read ()
{
   cout << "Reading input file " << input_file << endl;
   Reader fin(input_file);

   read_grid (fin);
   read_numeric (fin);
   read_material (fin);
   read_constants (fin);
   read_initial_condition (fin);
   read_exact (fin);
   read_boundary_condition (fin);
   read_output (fin);
}

//------------------------------------------------------------------------------
// Read grid section
//------------------------------------------------------------------------------
void Parameter::read_grid (Reader &fin)
{
   cout << "  Reading grid section\n";

   string input;
   int N;

   fin.begin_section ("grid");

   fin.entry ("x_min");
   fin.getline (input);
   domain.add ("xmin", input);
   fin.entry ("x_max");
   fin.getline (input);
   domain.add ("xmax", input);
   

   fin.begin_section ("Nx_cell");
   while(!fin.eos())
   {
      fin >> N;
      assert(N > 0);
      Nx_cell.push_back(N);
   }
   
   
   fin.entry ("y_min");
   fin.getline (input);
   domain.add ("ymin", input);
   fin.entry ("y_max");
   fin.getline (input);
   domain.add ("ymax", input);
   fin.begin_section ("Ny_cell");
   while(!fin.eos())
   {
      fin >> N;
      assert(N > 0);
      Ny_cell.push_back(N);
   }
   
   fin.entry ("z_min");
   fin.getline (input);
   domain.add ("zmin", input);
   fin.entry ("z_max");
   fin.getline (input);
   domain.add ("zmax", input);
   fin.begin_section ("Nz_cell");
   while(!fin.eos())
   {
      fin >> N;
      assert(N > 0);
      Nz_cell.push_back(N);
   }

   domain.value(xmin,xmax,ymin,ymax,zmin,zmax);
   assert(xmin<xmax);
   assert(ymin<ymax);
   assert(zmin<zmax);

   fin.end_section ();
   
   assert(Nx_cell.size() == Ny_cell.size());
   assert(Nx_cell.size() == Nz_cell.size());
   
}

//------------------------------------------------------------------------------
// Read numeric section
//------------------------------------------------------------------------------
void Parameter::read_numeric (Reader &fin)
{
   cout << "  Reading numeric section\n";

   string input;

   fin.begin_section ("numeric");

   fin.entry ("time_scheme");
   fin >> input;
   if(input=="rk1")
   {
      n_rks = 1;
      time_scheme = rk1;
   }   
   else if(input=="ssprk3")
   {
      n_rks = 3;
      time_scheme = ssprk3;
   } 
   else if(input=="jameson_rk4")
   {
      n_rks = 4;
      time_scheme = jameson_rk4;
   }
   else if(input=="rk4")
   {
      n_rks = 4;
      time_scheme = rk4;
   } 
   else
   {
      cout << "   Error: unknown time_scheme " << input << endl;
      exit (0);  
   }

   fin.entry ("cfl");
   fin >> cfl;
   assert (cfl > 0.0);

   fin.entry ("max_iter");
   fin >> max_iter;
   assert (max_iter > 0);

   fin.entry ("final_time");
   fin >> final_time;
   assert (final_time >= 0.0);

   fin.entry ("min_residue");
   fin >> min_residue;
   assert (min_residue > 0.0);

   fin.entry("reconstruct");
   fin >> input;
   reconstruction = input;
   if(input == "first")
   {
      reconstruct_scheme = first;
      rec_ng_cell = 1;
   }   
   else if(input == "minmod") 
   {  
      reconstruct_scheme = tvd_minmod;
      rec_ng_cell = 2;
   }
   else if(input == "minabs") 
   {  
      reconstruct_scheme = minabs_lim;
      rec_ng_cell = 2;
   }
   else if(input == "second") 
   {  
      reconstruct_scheme = second;
      rec_ng_cell = 2;
   } 
   else if(input == "eno") 
   {  
      reconstruct_scheme = eno;
      fin >> eno_level;
      rec_ng_cell = eno_level;
      reconstruction += to_string(eno_level);
   }
   else if(input == "sp_weno") 
   {  
      reconstruct_scheme = sp_weno;
      spweno_corr_type = 0;
      rec_ng_cell = 2;
   } 
   else if(input == "sp_wenoc") 
   {  
      reconstruct_scheme = sp_weno;
      spweno_corr_type = 1;
      rec_ng_cell = 2;
   } 
   else if(input == "sp_weno_corr2") 
   {  
      reconstruct_scheme = sp_weno;
      spweno_corr_type = 2;
      rec_ng_cell = 2;
   }
   else if(input == "sp_weno_dl")
   {
      reconstruct_scheme = sp_weno_dl;
      rec_ng_cell = 2;
   } 
   else
   {
      cout << "   Error: unknown reconstruction method " << input << endl;
      exit (0);   
   }
   fin.entry ("net_file");
   fin >> net_file;

   fin.end_section ();
}    

//------------------------------------------------------------------------------
// Read material section
//------------------------------------------------------------------------------
void Parameter::read_material (Reader &fin)
{
   cout << "  Reading material section\n";

   string input;

   fin.begin_section ("material");

   fin.entry ("gamma");
   fin >> GAMMA;
   assert (GAMMA > 1.0);

   fin.entry ("gas_const");
   fin >> gas_const;
   assert (gas_const > 0.0);

   fin.begin_section ("viscosity");
   fin.entry ("model");
   fin >> input;
   if(input == "constant")
   {
      mu_model = mu_constant;

      fin.entry ("mu_ref");
      fin >> mu_ref;
      assert (mu_ref >= 0.0);
   }
   else if(input == "sutherland")
   {
      mu_model = mu_sutherland;
      fin.entry ("mu_ref");
      fin >> mu_ref;
      fin.entry ("T_ref");
      fin >> T_ref;
      fin.entry ("T_0");
      fin >> T_0;
      assert (mu_ref > 0.0);
      assert (T_ref > 0.0);
      assert (T_0 > 0.0);
   }
   else if(input == "power")
   {
      mu_model = mu_power;
      fin.entry ("mu_ref");
      fin >> mu_ref;
      fin.entry ("T_ref");
      fin >> T_ref;
      fin.entry ("omega");
      fin >> omega;
      assert (mu_ref > 0.0);
      assert (T_ref > 0.0);
      assert (omega > 0.0);
   }
   else
   {
      cout << "   Error: unknown viscosity type " << input << endl;
      exit (0);
   }
   fin.end_section ();

   fin.entry ("prandtl");
   fin >> Pr;
   assert (Pr > 0.0);
   
   Cp = GAMMA * gas_const / (GAMMA - 1.0);

   fin.entry ("model");
   fin >> input;
   if(input == "euler")
      model = euler;
   else if(input == "ns_kep")
      model = ns_kep;
   else if(input == "ns_es")
      model = ns_es;   
   else
   {
      cout << "   Error: unknown flow model " << input << endl;
      exit (0);
   }
   
   if(mu_ref == 0.0 && (model == ns_es || model == ns_kep))
   {
      cout << "   Error: cannot set mu_ref = 0 for Navier-Stokes " << endl;
      exit (0);
   }

   fin.entry ("flux");
   fin >> input;
   if(input == "kepes_tecno_roe")
   {
      flux_scheme = kepes_tecno_roe;  
      flux_ng_cell = 1;
   }
   else if(input == "kepes_tecno4_roe")
   {
      flux_scheme = kepes_tecno4_roe;  
      flux_ng_cell = 2;
   }   
   else if(input == "kepes_tecno_rusanov")
   {
      flux_scheme = kepes_tecno_rusanov;
      flux_ng_cell = 1;   
   }
   else if(input == "kepes_tecno4_rusanov")
   {
      flux_scheme = kepes_tecno4_rusanov;
      flux_ng_cell = 2;   
   }   
   else if(input == "kepec")
   {
      flux_scheme = kepec;
      flux_ng_cell = 1;
   } 
   else if(input == "kepec4")
   {
      flux_scheme = kepec4;
      flux_ng_cell = 2;
   } 
   else if(input == "kep")
   {
      flux_scheme = kep;
      flux_ng_cell = 1;
   }
   else if(input == "kep4")
   {
      flux_scheme = kep4;
      flux_ng_cell = 2;
   }
   else if(input == "keps")
   {
      flux_scheme = keps;
      flux_ng_cell = 1;
   }
   else if(input == "keps4")
   {
      flux_scheme = keps4;
      flux_ng_cell = 2;
   }
   else if(input == "kep_tecno_roe")
   {
      flux_scheme = kep_tecno_roe;  
      flux_ng_cell = 1;
   }
   else if(input == "roe")
   {
      flux_scheme = roe;
      flux_ng_cell = 1;
   } 
   else if(input == "roe_ec")
   {
      flux_scheme = roe_ec;
      flux_ng_cell = 1;
   } 
   else if(input == "roe_ec4")
   {
      flux_scheme = roe_ec4;
      flux_ng_cell = 2;
   }  
   else if(input == "roe_tecno_roe")
   {
      flux_scheme = roe_tecno_roe;  
      flux_ng_cell = 1;
   }
   else if(input == "roe_tecno4_roe")
   {
      flux_scheme = roe_tecno4_roe;  
      flux_ng_cell = 2;
   } 
   else   
   {
      cout << "   Error:: unknown flux scheme: " << input << endl;
      exit (0);
   }
   
   fin.entry("EC1_fix");   
   fin >> input;
   if(input == "yes")
      EC1_fix = true;
   else if(input == "no")
      EC1_fix = false;   
   else   
   {
      cout << "   Error:: unknown EC1_fix option: " << input << endl;
      exit (0);
   }
   
   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read constants section
//------------------------------------------------------------------------------
void Parameter::read_constants (Reader &fin)
{
   cout << "  Reading constants section\n";

   string input;
   double value;

   fin.begin_section ("constants");

   while (!fin.eos())
   {
      fin >> input;
      fin >> value;
      cout << setw(16) << input << setw(16) << value << endl;
      constants.insert ( pair<string,double>(input, value) );
   }

}

//------------------------------------------------------------------------------
// Read initial condition
//------------------------------------------------------------------------------
void Parameter::read_initial_condition (Reader &fin)
{
   cout << "  Reading initial condition section\n";

   string input;

   fin.begin_section ("initial_condition");
      
   fin.entry("density");
   fin.getline (input);
   initial_condition.add ("density", input);

   fin.entry ("xvelocity");
   fin.getline (input);
   initial_condition.add ("xvelocity", input);
  
   fin.entry ("yvelocity");
   fin.getline (input);
   initial_condition.add ("yvelocity", input);
   
   fin.entry ("zvelocity");
   fin.getline (input);
   initial_condition.add ("zvelocity", input);

   fin.entry ("pressure");
   fin.getline (input);
   initial_condition.add ("pressure", input);

   
   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read exact solution
//------------------------------------------------------------------------------
void Parameter::read_exact (Reader &fin)
{
   cout << "  Reading exact solution section\n";

   string input;

   fin.begin_section ("exact_soln");
   
   fin.entry("available");   
   fin >> input;
   if(input == "yes")
      exact_available = true;
   else if(input == "no")
      exact_available = false;   
   else   
   {
      cout << "   Error:: unknown exact solution availability option: " << input << endl;
      exit (0);
   }
   
   if(exact_available)
   {   
	  fin.entry("density");
	  fin.getline (input);
	  exact_soln.add ("density", input);

	  fin.entry ("xvelocity");
	  fin.getline (input);
	  exact_soln.add ("xvelocity", input);
  
	  fin.entry ("yvelocity");
	  fin.getline (input);
	  exact_soln.add ("yvelocity", input);
	  
	  fin.entry ("zvelocity");
	  fin.getline (input);
	  exact_soln.add ("zvelocity", input);

	  fin.entry ("pressure");
	  fin.getline (input);
	  exact_soln.add ("pressure", input);
	  
	  fin.end_section ();
   }
   else
   {
      while(!fin.eos())
        fin >> input;
   }     
}

//------------------------------------------------------------------------------
// Read boundary conditions
//------------------------------------------------------------------------------
void Parameter::read_boundary_condition (Reader &fin)
{
   cout << "  Reading boundary condition section\n";
   string input;
   fin.begin_section ("boundary_condition");

   string boundary_names[] = {"xleft","xright","yleft","yright","zleft","zright"};
   for(int i=0; i<6;++i)
   {
      fin.begin_section (boundary_names[i]);
      fin.entry("type");
      fin >> input;
	  if(input=="periodic")
	  {
		 bc[i] = periodic;
	  } 
	  else if(input=="dirichlet")
	  {
		 bc[i] = dirichlet;
	  } 
	  else if(input=="neumann")
	  {
		 bc[i] = neumann;
	  } 
	  else if(input=="wall")
	  {
		 bc[i] = wall;
	  } 
	  else
	  {
		 cout << "   Error: unknown boundary type " << input << endl;
		 exit (0);  
	  } 
	  if(input == "dirichlet")
	  {
		 cout << "   Error: Dirichlet type BC not implemented yet "<<endl;
		 exit (0);  
	  }
      fin.end_section();
      if((model==ns_kep || model==ns_es) && input!="periodic")
      {
		 cout << "   Error: Only periodic BC allowed with Navier-Stokes for the time being "<<endl;
		 exit (0);  
	  }
   }
   if(bc[0] == periodic)
     assert(bc[1] == periodic);
   if(bc[1] == periodic)
     assert(bc[0] == periodic);
   if(bc[2] == periodic)
     assert(bc[3] == periodic);   
   if(bc[3] == periodic)
     assert(bc[2] == periodic);
   if(bc[4] == periodic)
     assert(bc[5] == periodic);   
   if(bc[5] == periodic)
     assert(bc[4] == periodic);     
   
   fin.end_section ();
   
   
}

//------------------------------------------------------------------------------
// Read output section
//------------------------------------------------------------------------------
void Parameter::read_output (Reader &fin)
{
   cout << "  Reading output section\n";

   string input;

   fin.begin_section ("output");
   
   fin.entry("OOC_study");   
   fin >> input;
   if(input == "yes")
      OOC = true;
   else if(input == "no")
      OOC = false;   
   else   
   {
      cout << "   Error:: unknown OOC_study option: " << input << endl;
      exit (0);
   }
   
   fin.entry("problem_dim");
   fin >> pdim;
   
   assert(pdim == 1 || pdim == 2 || pdim == 3);
   
   if(OOC && (Nx_cell.size() <2 || Ny_cell.size() < 2))
   {
      cout << "   Error:: Cannot perform OOC study without at least two mesh levels"<<endl;
      exit (0);
   }
   
   if(OOC && !exact_available)
   {
      cout << "   Error:: Cannot perform OOC study as exact solution not avaible"<<endl;
      exit (0);
   }
   
   fin.entry ("time_stamps");
   fin >> t_stamps;
   assert (t_stamps >= 1);
   
   fin.entry ("output_path");
   fin >> OUTPUT_PATH;

   fin.end_section ();
   
   if(OOC)
      mesh_levels = Nx_cell.size();
   else
      mesh_levels = 1;   
}

//------------------------------------------------------------------------------
// Viscosity coefficient
//------------------------------------------------------------------------------
double Parameter::Viscosity (const double T) const
{
   switch (mu_model)
   {
      case mu_constant:
         return mu_ref;

      case mu_sutherland:
         return mu_ref * pow(T/T_ref, 1.5) * (T_ref + T_0) / (T + T_0);

      case mu_power:
         return mu_ref * pow(T/T_ref, omega);
         
      default:
         cout << "    viscosity: unknown model " << mu_model << endl;
         abort ();
   }
}
