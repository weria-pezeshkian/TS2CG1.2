#if !defined(AFX_SurfaceMap_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_SurfaceMap_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

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
#include "molecules.h"
#include "point.h"
#include "inclusion.h"
#include "GenerateMolType.h"
#include "Tensor2.h"
#include "ReadDTSFolder.h"

struct ProteinList {
    std::string ProteinName;
    double Ratio;
    double Phi;
    double Theta;
    double Z0;
    int ID;
    int Maxno;
    int created;
    
} ;
struct ExcludedVolumeBeads {
    Vec3D X;
    double R;
    
} ;
class SurfaceMap
{
public:
    
	SurfaceMap(Argument *pArgu);
	virtual ~SurfaceMap();
    



private:
    Vec3D  m_Box;
    std::vector<point>  m_PointUp;
    std::vector<point>  m_PointDown;
    std::vector<point>  m_WallPointUp;
    std::vector<point>  m_WallPointDown;
    std::vector<inclusion>  m_Inclusion;

    
public:
    Vec3D  *m_pBox;
    std::vector<point*>  m_pPointUp;
    std::vector<point*>  m_pPointDown;
    std::vector<point*>  m_pWallPointUp;
    std::vector<point*>  m_pWallPointDown;
    std::vector<inclusion*>  m_pInclusion;
    
    
    
    std::vector<inclusion*>  m_pInc;
    std::string m_InclusionDirectionType;
    std::vector<exclusion*>  m_pExc;

private:
    double dist2between2Points(Vec3D X1,Vec3D X2);
    bool AnyBeadWithinCutoff(UnitCell *,Vec3D Pos);
    Tensor2 Rz(double cos, double sin);
    Tensor2 TransferMatLG(Vec3D N, Vec3D t1, Vec3D t2);
    double    m_AvailPointDown ;
    bool m_monolayer;
    std::vector<point*>  m_point1;
    std::vector<point*>  m_point2;


    std::string functiontype(std::string filename);
};


#endif
