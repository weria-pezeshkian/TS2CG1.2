#if !defined(AFX_Membrane_Builder_H_884B21B8_C13D_5648_BF23_124095086234__INCLUDED_)
#define AFX_Membrane_Builder_H_884B21B8_C13D_5648_BF23_124095086234__INCLUDED_

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
#include "SHGeneric1DPBCPointMaker.h"

class Membrane_Builder
{
public:
    
	Membrane_Builder(Argument *pArgu);
	virtual ~Membrane_Builder();
    


public:

private:
    std::vector<bead*> m_pAllBeads;
    std::vector<molecules> m_AllMolecules;
    std::vector<molecules*> m_pAllMolecules;
    Vec3D *m_pBox;
    Vec3D  m_Box;
    std::vector<bead> m_FinalBeads;
    std::string m_FinalOutputGroFileName;
    std::string m_FinalTopologyFileName;

    std::map<std::string , MolType>  m_MoleculesType;
    int m_ResID;
    double m_Iter;

private:
    double dist2between2Points(Vec3D X1,Vec3D X2);
    void WriteFinalGroFile();
    Tensor2 Rz(double cos, double sin);
    Tensor2 TransferMatLG(Vec3D N, Vec3D t1, Vec3D t2);
    bool m_Renormalizedlipidratio;
    double    m_TotalAreaUp ;
    double    m_TotalAreaDown ;
    double    m_AvailPointUp ;
    double    m_AvailPointDown ;
    double    m_APLLipids ;
    double    m_APLWall ;


    bool m_monolayer;
    std::vector<point>  m_Point1;
    std::vector<point>  m_Point2;
    std::vector<point>  m_WallPoint1;
    std::vector<point>  m_WallPoint2;
    std::vector<point*>  m_pWallPoint1;
    std::vector<point*>  m_pWallPoint2;
    std::vector<point*>  m_pPoint1;
    std::vector<point*>  m_pPoint2;
   // void CreateWallBead(std::vector<point*>  p1, std::vector<point*>  p2);
    void GenLipid(MolType moltype, int , Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2);
    std::string functiontype(std::string filename); /// Read the data from the str file.





};


#endif
