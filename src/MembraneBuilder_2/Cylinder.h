#if !defined(AFX_Cylinder_H_884B21B8_C13D_5648_BF23_124095086234__INCLUDED_)
#define AFX_Cylinder_H_884B21B8_C13D_5648_BF23_124095086234__INCLUDED_

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
#include "point.h"
#include "inclusion.h"
#include "GenerateMolType.h"
#include "Tensor2.h"


class Cylinder
{
public:
    
	Cylinder(Argument *pArgu);
	virtual ~Cylinder();
    
    inline  std::vector<point> GetWallPoint1()                const  {return m_WallPoint1;} // returns all the created wall beads
    inline  std::vector<point> GetWallPoint2()                const  {return m_WallPoint2;} // returns all the created wall beads

    inline  std::vector<point> GetUpPoint()                const  {return m_Point1;}    // upper monolayer beads
    inline  std::vector<point> GetInPoint()                const  {return m_Point2;}    // inner monolayer beads

    inline  Vec3D GetBox()                const  {return m_Box;}    // returns the box sides


public:

private:
    Vec3D *m_pBox;
    Vec3D  m_Box;   // system box size
    double m_Thickness;     // thickness of the membrane
    double m_Density;       // density of the lipids points, it should be smaller then the smaller area per lipids of the select lipids
    double m_WallDensityup;   // density of the wall beads. nm^-2
    double m_WallDensityin;   // density of the wall beads. nm^-2
    double m_DL;   // density of the wall beads. nm^-2
    double m_R;   // density of the wall beads. nm^-2

    std::vector<double>  m_WallBox;         // a volume to have the wall beads in.

private:
    void Initialize(std::string filename); /// Read the data from the str file.

    std::vector<point> CalculateArea_MakePoints(int i, double APL,double H,double DL, bool); // a function to create the area and the beads


    bool m_monolayer;
    std::vector<point>  m_Point1;
    std::vector<point>  m_Point2;
    std::vector<point>  m_WallPoint1;
    std::vector<point>  m_WallPoint2;
    Vec3D F(double t,int layer, double H);
    Vec3D Normal(double t);



};


#endif
