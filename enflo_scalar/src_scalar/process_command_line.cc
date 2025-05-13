#include <iostream>
#include <cassert>
#include <cstring>
#include <cstdlib>

using namespace std;

void show_options ();

//------------------------------------------------------------------------------
// Get command line flags and input file
//------------------------------------------------------------------------------
void process_command_line (int   argc,
                           char* argv[],
                           int&  ifile)
{ 
   if(argc < 2)
      show_options ();

   int i = 1;
   bool found_input_file = false;

   while (i < argc)
   {
      
      if(strcmp(argv[i],"-i")==0)
      {
         assert (i+1 < argc); 
         ifile = i+1;
         ++i;
         found_input_file = true;
      }
      else
      {
         cout << "Unknown command line flag: " << argv[i] << endl;
         show_options ();
      }

      ++i;
   }

   if(!found_input_file)
      show_options ();
}

//------------------------------------------------------------------------------
// Print command line options available
//------------------------------------------------------------------------------
void show_options ()
{

   cout << "Valid flags are:\n";
   cout << "   -i filename   Specify input file name (required)\n";

}
