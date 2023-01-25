#include <iostream>
#include "help.h"
#include "Def.h"

help::help(std::string exe)
{
    int size= exe.size();
 {
     

    std::cout<<"\n";
     std::cout<<"=========================================================================================================="<<"\n";
     std::cout<<"-- TS2CG solvate (SOL)   "<<"\n";
     std::cout<<"-- Version:  "<<SoftWareVersion<<"\n";
     std::cout<<"-- Groningen Biomolecular Sciences and Biotechnology Institute and Zernike Institute for Advanced Materials,\n-- University of Groningen, Groningen, Netherlands"<<"\n";
     std::cout<<"-- For more information contact Weria Pezeshkian: w.pezeshkian@rug.nl, weria.pezeshkian@gmail.com"<<"\n";
     std::cout<<"-- citation: Pezeshkian, et al. Nat. Comm. 11, 2296 (2020)."<<"\n";
     std::cout<<"=========================================================================================================="<<"\n";
     std::cout<<"-- With option -Bondlength, you can chnage the initial bond guess. Large Bondlength may generate an unstable structure ";

     std::cout<<"=========================================================================================================="<<"\n";
    std::cout<<"------------ This script convert Pointillism outputs to a CG model -------------------"<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  option            type        default            description "<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  -in              string       input.gro          name of the file to be solvated "<<"\n";
    std::cout<<"  -o               string       output.gro         output file name "<<"\n";
    std::cout<<"  -Rcutoff          double       0.4                cutoff distance "<<"\n";
    std::cout<<"  -tem             string       W.gro              name of the template file for solvation"<<"\n";


    

     std::cout<< "example: SOL  "<<"\n";




   }

    

}

help::~help()
{
    
}
