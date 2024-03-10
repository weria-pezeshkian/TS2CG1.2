#include <iostream>
#include "help.h"
help::help(int version,std::string exe)
{
    int size= exe.size();
 {

    std::cout<<"\n";
     std::cout<<"------------This software was developed by Weria Pezeshkian at Univeristy of Groningen \n Copyright  Weria Pezeshkian -------------------"<<"\n";
     std::cout<<"------------ any question: contact weria.pezeshkian@gmail.com -------------------"<<"\n";
    std::cout<<"------------simple example for exacuting clustering analysis -------------------"<<"\n";
    std::cout<< "./C++Analysis -xtc file.xtc  -b 0  -e 100 -gro file.gro  -ndx file.ndx -cutoff 1.1  "<<"\n";
     std::cout<<"------------ and for flipflop run -------------------"<<"\n";
     std::cout<< "./C++Analysis -xtc file.xtc  -b 0  -e 100 -gro file.gro  -ndx file.ndx -ndx2 file2.ndx -cutoff 1.1 -gangle 20 -r flipflop "<<"\n";
     std::cout<<"------------ and for order parameter run -------------------"<<"\n";
     std::cout<< "./C++Analysis -xtc file.xtc  -b 0  -e 100 -gro file.gro  -ndx file.ndx  -r orderparameter -bins 200 -tf 100 "<<"\n";




    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  option    type        default            description "<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  -xtc       string       MD.xtc                input file name "<<"\n";
    std::cout<<"  -b         int          0      	            initial time step "<<"\n";
    std::cout<<"  -e         int          max	                final time step  "<<"\n";
    std::cout<<"  -r         string       clustering            analsysis type (clustering/flipflop)  "<<"\n";
    std::cout<<"  -seed      int          737      	            Random number seed  "<<"\n";
    std::cout<<"  -defout    string       output      	        an string for run out put files   "<<"\n";
    std::cout<<"  -ndx       string       Index.ndx             Index file name  "<<"\n";
    std::cout<<"  -ndx2       string      Index2.ndx            second Index file name  "<<"\n";
    std::cout<<"  -cutoff     double      1.1                   cutoff  "<<"\n";
    std::cout<<"  -gangle     double      0                     an angle to define middle of the bilayer  "<<"\n";
    std::cout<<"  -more       string      no(yes)               print in files more data  "<<"\n";
    std::cout<<"  -tf         int         10                    no frame to make average of  "<<"\n";
     std::cout<<"  -rg         string     DOPC                  reference group  "<<"\n";
     std::cout<<"  -tg         int        PO4                   target group group  "<<"\n";
     std::cout<<"  -bl         double     0.1                   bin size  "<<"\n";

    std::cout<<"  -bins       int         10                     divide the box into X bins "<<"\n";

    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"=========================================================================="<<"\n";
     std::cout<<"------------------ version: "<< version<< " ------------------"<<"\n";
    std::cout<<" =================================================================  \n";


   }

    

}

help::~help()
{
    
}
