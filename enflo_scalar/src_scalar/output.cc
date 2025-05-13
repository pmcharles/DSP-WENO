#include <algorithm>
#include <cmath>
#include "fv.h"

using namespace std;

//------------------------------------------------------------------------------
// Solution output
//------------------------------------------------------------------------------
void FV::output ()
{
   
   string filename = "sol";
   string filename2 = "1Dsol";
   string extension = ".vtk";
   string extension2 = ".dat";

   string precount;
   if     (counter <= 9)    precount = "000";
   else if(counter <= 99)   precount = "00";
   else if(counter <= 999)  precount = "0";
   else if(counter <= 9999) precount = "";
   else
   {
      cout << "Writer::output: counter is too large !!!\n";
   }
   
   stringstream ss;
   ss <<precount<<counter;
   filename += ss.str();
   filename +=extension;
   filename = SOLN_DIR + "/" + filename;
   filename2 += ss.str();
   filename2 +=extension2;
   filename2 = SOLN_DIR + "/" + filename2;
    
   cout<<"Saving solutions in "<< filename <<endl;   
 
   ofstream vtk(filename.c_str());
   
   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "enflo" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET STRUCTURED_POINTS" << endl;
   vtk << "FIELD FieldData 1" <<endl;
   vtk << "TIME 1 1 double"<<endl;
   vtk << time <<endl;
   vtk << "DIMENSIONS " << (Nx_cell+1) << " " << (Ny_cell+1) << " " << (Nz_cell+1) << endl;
   vtk << "SPACING " << dx << " " << dy << " " << dz << " "  << endl;
   vtk << "ORIGIN " << param.xmin << " " << param.ymin << " " << param.zmin << " "  << endl;

   vtk << "CELL_DATA  " << (Nx_cell)*(Ny_cell)*(Nz_cell) << endl;

   vtk << "SCALARS solution double 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(int k=0; k<Nz_cell; ++k)
      for(int j=0; j<Ny_cell; ++j)
         for(int i=0; i<Nx_cell; ++i)
             vtk << solution[i+ng_cell][j+ng_cell][k+ng_cell] << endl;         
                       
   vtk.close ();
   
   ofstream file(filename2.c_str());
   for (int i=0; i<Nx_cell; ++i)
      file << (param.xmin+dx/2)+dx*i << "\t" << solution[i+ng_cell][ng_cell][ng_cell] << "\t" << endl;

   counter++;
}
