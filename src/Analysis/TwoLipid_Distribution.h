#if !defined(AFX_TwoLipid_Distribution_H_7F4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_TwoLipid_Distribution_H_7F4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

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

class TwoLipid_Distribution
{
public:
    
	TwoLipid_Distribution(Argument *pArgu);
	virtual ~TwoLipid_Distribution();
    


public:

private:
    std::vector<bead*> m_pAllBeads;
    std::vector<bead*> m_pWallBeads;
    std::vector<molecules> m_AllMolecules;
    std::vector<molecules*> m_pAllMolecules;
    Vec3D *m_pBox;
    std::vector<double> FindCurv (bead *p);
    
    
private:
    double dist2between2bead(bead*,bead*);
    Vec3D COG_PBC(Vec3D *pBox, std::vector<bead *> pBeads);
    std::vector<point> m_WallPoints;

};


#endif
