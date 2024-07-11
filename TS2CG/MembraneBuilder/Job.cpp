


#include "Argument.h"
#include "Job.h"
#include "Nfunction.h"
#include "BackMap.h"
// this class does not do much, it is just an extra check in case in future we want to diversify.
// this class just get the arguments and call BackMap class. 
Job::Job(std::vector<std::string> argument) {
    Argument arg(argument);
    std::string executable = SubstringFromRight(argument.at(0));
    std::string function = arg.GetFunction();
//-- checking if the ExecutableName is PCG and also if the type of fuctions known,
    if (executable==ExecutableName) { // ExecutableName = PCG
        if (function == "backmap" || function == "analytical_shape") {
            BackMap B(&arg); // call Backmap class
        } else {
            std::cout << function << "---> function is not recognized \n";
            std::exit(0);
        }
    }
    else {
        std::cout << "---> error executable name <"<<executable<<">  is wrong \n";
        std::exit(0);
    }
}
Job::~Job()
{
    
}
std::string Job::SubstringFromRight(const std::string& input) {
    size_t pos = input.find_last_of('/');
    if (pos != std::string::npos) {
        return input.substr(pos + 1);
    } else {
        // If there is no '/' in the string, return the entire string
        return input;
    }
}
