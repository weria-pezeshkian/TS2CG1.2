/* This code was developed by Weria Pezeshkian at the University of Groningen
 Copyright (c) Weria Pezeshkian
 */
#include "Def.h"
#include "Job.h"
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char* argv[])
{
    std::vector<std::string> argument;
    std::string Temprory;

    if (argc == 1)
    {
        // This means, the user wants to run the command without any argument
        std::cout << "Error: not enough arguments to perform the analysis";
        return 0;
    }
    else
    {
        for (long i = 0; i < argc; i++)
        {
            Temprory.assign(argv[i]);
            argument.push_back(Temprory);
        }

        Job job(argument);
    }

    return 0;
}
