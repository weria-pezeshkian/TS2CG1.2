#include <iostream>
#include "help.h"
#include "Def.h"

help::help(std::string exe) {

    std::cout << "\n";
    std::cout << "==========================================================================================================" << "\n";
    std::cout << "--> " << SoftWareName << "\n";
    std::cout << "-- Version:  " << SoftWareVersion << "\n";
    std::cout << "-- Niels Bohr Institute,  \n University of Copenhagen, Copenhagen, Denmark" << "\n";
    std::cout << "-- For more information contact Weria Pezeshkian: weria.pezeshkian@nbi.ku.dk, weria.pezeshkian@gmail.com" << "\n";
    std::cout << "--> citation: Pezeshkian, et al. Nat. Comm. 11, 2296 (2020)." << "\n";
    std::cout << "==========================================================================================================" << "\n";

    std::cout << "==========================================================================================================" << "\n";
    std::cout << "------------ This script solvates the system -------------------" << "\n";
    std::cout << "-------------------------------------------------------------------------------" << "\n";
    std::cout << "  option            type        default            description " << "\n";
    std::cout << "-------------------------------------------------------------------------------" << "\n";
    std::cout << "  -in              string       input.gro          name of the file to be solvated " << "\n";
    std::cout << "  -ion             integers     0 0                generating ions; two integers number of positive and negativer ions" << "\n";
    std::cout << "  -db              double       0.05               the distance between each copy when propagating the water box" << "\n";
    std::cout << "  -o               string       output.gro         output file name " << "\n";
    std::cout << "  -Rcutoff         double       0.4                cutoff distance " << "\n";
    std::cout << "  -tem             string       W.gro              name of the template file for solvation" << "\n";
    std::cout << "  -nname           string      CL                  name of the negative ions" << "\n";
    std::cout << "  -pname           string      NA                  name of the positive ions" << "\n";
    std::cout << "  -seed            integer     9474                seed for ion placement " << "\n";
    std::cout << "  -unsize          double      2                   size of the unitcells for overlap checking, smaller numbers are faster but needs more RAM, should not be smaller than Rcutoff  " << "\n";
    std::cout << "example: " << "\n";
    std::cout << ExcName << "   -in in.gro -o out.gro -ion 20 20  -tem water.gro" << "\n";
}

help::~help() {

}
