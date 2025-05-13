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

   vtk << "SCALARS density float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(int k=0; k<Nz_cell; ++k)
      for(int j=0; j<Ny_cell; ++j)
         for(int i=0; i<Nx_cell; ++i)
             vtk << primitive[i+ng_cell][j+ng_cell][k+ng_cell].density << endl;
             
   vtk << "VECTORS velocity float" << endl;
   for(int k=0; k<Nz_cell; ++k)
      for(int j=0; j<Ny_cell; ++j)
         for(int i=0; i<Nx_cell; ++i)
             vtk << primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity.x << " "
                 << primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity.y << " "
                 << primitive[i+ng_cell][j+ng_cell][k+ng_cell].velocity.z << endl;          
   
   vtk << "SCALARS pressure float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(int k=0; k<Nz_cell; ++k)
      for(int j=0; j<Ny_cell; ++j)
         for(int i=0; i<Nx_cell; ++i)
             vtk << primitive[i+ng_cell][j+ng_cell][k+ng_cell].pressure << endl;
             
   vtk << "SCALARS physical_entropy float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(int k=0; k<Nz_cell; ++k)
      for(int j=0; j<Ny_cell; ++j)
         for(int i=0; i<Nx_cell; ++i)
             vtk << Entropy(primitive[i+ng_cell][j+ng_cell][k+ng_cell]) << endl; 
             
   
   
  //  PrimVar exact_pw;
//    
//    vtk << "SCALARS Edensity float 1" << endl;
//    vtk << "LOOKUP_TABLE default" << endl;
//    for(int k=0; k<Nz_cell; ++k)
//       for(int j=0; j<Ny_cell; ++j)
//          for(int i=0; i<Nx_cell; ++i)
//          {   
//              param.exact_soln.value (cell_cc[i][j][k],time,exact_pw);
//              vtk << exact_pw.density << endl;
//          }    
//              
//    vtk << "VECTORS Evelocity float" << endl;
//    for(int k=0; k<Nz_cell; ++k)
//       for(int j=0; j<Ny_cell; ++j)
//          for(int i=0; i<Nx_cell; ++i)
//          {
//              param.exact_soln.value (cell_cc[i][j][k],time,exact_pw);
//              vtk << exact_pw.velocity.x << " "
//                  << exact_pw.velocity.y << " "
//                  << exact_pw.velocity.z << endl; 
//          }                 
//    
//    vtk << "SCALARS Epressure float 1" << endl;
//    vtk << "LOOKUP_TABLE default" << endl;
//    for(int k=0; k<Nz_cell; ++k)
//       for(int j=0; j<Ny_cell; ++j)
//          for(int i=0; i<Nx_cell; ++i)
//          {
//              param.exact_soln.value (cell_cc[i][j][k],time,exact_pw);
//              vtk << exact_pw.pressure << endl;
//          }    
//              
//    vtk << "SCALARS Ephysical_entropy float 1" << endl;
//    vtk << "LOOKUP_TABLE default" << endl;
//    for(int k=0; k<Nz_cell; ++k)
//       for(int j=0; j<Ny_cell; ++j)
//          for(int i=0; i<Nx_cell; ++i)
//          {
//              param.exact_soln.value (cell_cc[i][j][k],time,exact_pw);
//              vtk << Entropy(exact_pw) << endl;   
//          }             
                       
   vtk.close ();
   
   
   filename = "omg";
   extension = ".vtk";
   
   filename += ss.str();
   filename +=extension;
   filename = SOLN_DIR + "/" + filename;
    
   cout<<"Saving vorticty in "<< filename <<endl;   
 
   ofstream vtk1(filename.c_str());
   
   vtk1 << "# vtk DataFile Version 3.0" << endl;
   vtk1 << "enflo" << endl;
   vtk1 << "ASCII" << endl;
   vtk1 << "DATASET STRUCTURED_POINTS" << endl;
   vtk1 << "FIELD FieldData 1" <<endl;
   vtk1 << "TIME 1 1 double"<<endl;
   vtk1 << time <<endl;
   vtk1 << "DIMENSIONS " << (Nx_cell+1) << " " << (Ny_cell+1) << " " << (Nz_cell+1) << endl;
   vtk1 << "SPACING " << dx << " " << dy << " " << dz << " "  << endl;
   vtk1 << "ORIGIN " << param.xmin << " " << param.ymin << " " << param.zmin << " "  << endl;

   vtk1 << "POINT_DATA  " << (Nx_cell+1)*(Ny_cell+1)*(Nz_cell+1) << endl;
             
   vtk1 << "VECTORS vorticity float" << endl;
   for(int k=0; k<Nz_cell+1; ++k)
      for(int j=0; j<Ny_cell+1; ++j)
         for(int i=0; i<Nx_cell+1; ++i)
         {
             vtk1 << vorticity[i][j][k].x << " "
                  << vorticity[i][j][k].y << " "
                  << vorticity[i][j][k].z << endl;       
         }   
                                
   vtk1.close ();  
   
   ofstream file(filename2.c_str());
   for (int i=0; i<Nx_cell; ++i)
      file << (param.xmin+dx/2)+dx*i << "\t" << primitive[i+ng_cell][ng_cell][ng_cell].density << "\t" << primitive[i+ng_cell][ng_cell][ng_cell].velocity.x << "\t" << primitive[i+ng_cell][ng_cell][ng_cell].pressure << endl;

   counter++;
}
