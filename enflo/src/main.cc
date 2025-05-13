#include <iostream>
#include "fv.h"


using namespace std;

map<string,double> constants;

void process_command_line (int argc, char* argv[], int& ifile);


int main (int argc, char* argv[])
{
   cout << "=================================================================================\n";
   cout << "   Starting 3D Compressible Euler solver on Cartesian mesh\n";   
   cout << "   ---  Author: Deep Ray \n";
   cout << "   ---  Date  : 5 April , 2016 \n";
   cout << "=================================================================================\n";
   
   int ifile;
   process_command_line (argc, argv, ifile);
   
   FV fv(argv[ifile]);
   fv.run();
   
   cout << "\n\n---------------Solver has finished!!----------------------\n"; 
   
   return 0;
}
