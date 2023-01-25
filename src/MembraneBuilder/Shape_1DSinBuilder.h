#if !defined(AFX_Shape_1DSinBuilder_H_884B21B8_C13D_5648_BF23_124095086234__INCLUDED_)
#define AFX_Shape_1DSinBuilder_H_884B21B8_C13D_5648_BF23_124095086234__INCLUDED_

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

/*struct ProteinList {
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
    
} ;*/
class Shape_1DSinBuilder
{
public:
    
	Shape_1DSinBuilder(Argument *pArgu);
	virtual ~Shape_1DSinBuilder();
    


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
    std::vector<point*>  m_pPoint1;
    std::vector<point*>  m_pPoint2;
    Shape_1DSin m_State;
   // void CreateWallBead(std::vector<point*>  p1, std::vector<point*>  p2);
    void GenLipid(MolType moltype, int , Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2);
    double CalculateArea_MakePoints(int i, double APL);

    Vec3D F(double t,int layer);
    Vec3D Normal(double t,int layer);
    Vec3D T1(double t,int layer);



};


#endif
