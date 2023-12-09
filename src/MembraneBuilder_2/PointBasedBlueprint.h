#if !defined(AFX_PointBasedBlueprint_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_PointBasedBlueprint_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

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
#include "point.h"
#include "inclusion.h"
#include "Tensor2.h"
#include "ReadDTSFolder.h"


class PointBasedBlueprint
{
public:
    
	PointBasedBlueprint(Argument *pArgu);
	virtual ~PointBasedBlueprint();
    


public:
    Wall *m_pWall;
    std::vector<point*>  m_pPointUp;
    std::vector<point*>  m_pPointDown;
    std::vector<exclusion*>  m_pExc;
    std::vector<inclusion*>  m_pInc;
    Vec3D *m_pBox;
    bool m_monolayer;

    
    
private:

    std::vector<point>  m_PointUp;
    std::vector<point>  m_PointDown;
    std::vector<point>  m_WPointUp;
    std::vector<point>  m_WPointDown;
    std::vector<inclusion>  m_Inc;
    std::vector<exclusion>  m_Exc;
    std::vector<point*>  m_pWPointUp;
    std::vector<point*>  m_pWPointDown;
    Vec3D  m_Box;
    Wall m_Wall;

private:
    std::string functiontype(std::string filename);
};


#endif
