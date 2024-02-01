#if !defined(AFX_GenerateMolType_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_GenerateMolType_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

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
#include "LMatrix.h"
#include "Nfunction.h"
#include "Argument.h"
#include "Vec3D.h"
#include "bead.h"
#include "ReadLipidLibrary.h"



class GenerateMolType
{
public:
    
	GenerateMolType(Argument *pArgu);
	~GenerateMolType();
    std::map<std::string , MolType>  GetMolType()        {return m_MoleculesType;}


public:
    void  Mol_Def_DPhospholipid(std::string, std::string BeadDEf, double APL);

private:
    std::map<std::string , MolType>  m_MoleculesType;
    double m_BondLenght;

private:
    void LiBMol();
    double MolAreaCal(MolType);  // this is only valid for gro files.
};


#endif
