#if !defined(AFX_Density_1D_H_7F4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_Density_1D_H_7F4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

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
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "Nfunction.h"
#include "Argument.h"
#include "Vec3D.h"
#include "bead.h"
#include "molecules.h"
#include "point.h"

class Density_1D
{
public:
    
	Density_1D(Argument *pArgu);
	virtual ~Density_1D();
    


public:

private:
    std::vector<bead*> m_pAllBeads;
    std::vector<bead*> m_pWallBeads;
    std::vector<bead*> m_pUpWallBeads;
    std::vector<bead*> m_pInWallBeads;
    std::vector<double> m_UpBeadCurve;
    std::vector<double> m_LowerBeadCurve;
    std::vector<double> m_SUpBeadCurve;
    std::vector<double> m_SLowerBeadCurve;
    std::vector<double> m_NUp;;
    std::vector<double> m_NDown;
    
    
    std::vector<molecules> m_AllMolecules;
    std::vector<molecules*> m_pAllMolecules;
    Vec3D *m_pBox;
    std::vector<double> FindCurv (bead *p);
    double m_DX;
    
private:

     LMatrix SplitDenisty(LMatrix a);

};


#endif
