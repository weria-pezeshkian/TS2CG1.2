#if !defined(AFX_Data_Structure_H)
#define AFX_Data_Structure_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <list>
#include <map>
#include <iomanip>
#include <valarray>
#include "bead.h"

struct MolType {
    std::vector<bead> Beads;
    std::string MolName;
    int beadnumber;
    double molarea;
} ;
struct DomainLipid {
    std::string Name;           // name of the lipid
    double Ap;                  // AP of this lipid structure
    double Ratio;               // how many/ratio in the domain
    int no_created;          // how many has been created
    int MaxNo;                  // how many of this type we should create
    
} ;
struct ProteinList {
    int ID; // some how matches the name in the str file and the id in the tsi file
    std::string ProteinName;
    double Ratio;
    double Phi;
    double Theta;
    double Z0;
    
    int Maxno;
    int created;
    
} ;
struct ExcludedVolumeBeads {
    Vec3D X;
    double R;
    
} ;
#endif
