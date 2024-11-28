/* This code was developed by Weria Pezeshkian at Univeristy of Southern Denmark
 
 Copyright (c) Weria Pezeshkian
 */
#include "SimDef.h"
#include "Job.h"




int main(int argc, char* argv[])
{
/// This part get the name of the input file as command line
//// for later use we set all of them in a vector
    std::vector <std::string> argument;
    std::string Temprory;

           for (long i=0;i<argc;i++)
           {
               Temprory.assign(argv[i]);
               argument.push_back(Temprory);
           }
  
         Job job(argument);


    return 0;

}
