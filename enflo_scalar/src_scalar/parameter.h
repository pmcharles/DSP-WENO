#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "reader.h"
#include "ext_constants.h"
#include "ic.h"
#include "exact.h"
#include "domain.h"

class Parameter
{
   public:
      char* input_file;
      std::string ic_file;
      std::string source_file;
      std::string net_file;
      
      bool ic_file_read;

      int n_rks;
      int max_iter;
      double cfl;
      double final_time;
      double min_residue;
      double xmin,xmax,ymin,ymax,zmin,zmax;
      std::vector<int> Nx_cell, Ny_cell, Nz_cell;
      int flux_ng_cell;
      int rec_ng_cell;
      
      BCType bc[6];
      
      Vector lin_vel;
      
      Domain domain;
      FlowModel model;
      FluxScheme flux_scheme;
      ReconstructionScheme reconstruct_scheme;
      int eno_level;
      TimeIntegrationScheme time_scheme;
      InitialCondition initial_condition;
      Exact exact_soln;
      std::string OUTPUT_PATH;
      std::string reconstruction;

      bool exact_available, OOC;
      int spweno_corr_type;
      int pdim;
      int mesh_levels;
      int t_stamps;
      void read ();

   private:
      void read_grid (Reader&);
      void read_constants (Reader&);
      void read_numeric (Reader&);
      void read_material (Reader&);
      void read_initial_condition (Reader&);
      void read_exact (Reader&);
      void read_boundary_condition (Reader&);
      void read_output (Reader&);
};

#endif
