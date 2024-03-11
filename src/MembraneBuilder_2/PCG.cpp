/* This code was developed by Weria Pezeshkian at Univeristy of Groningen,
 updated in 2023 at the Univeristy of Copenhagen
 Copyright (c) Weria Pezeshkian
 email: weria.pezeshkian@gmail.com
 */
#include "Def.h"
#include "Job.h"

int main(int argc, char* argv[])
{
    // getting the commandline strings and passing them to job class
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
