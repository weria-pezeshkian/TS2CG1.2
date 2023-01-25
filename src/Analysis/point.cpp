

#include <stdio.h>
#include "point.h"
//#include "links.h"
point::point(int id,Vec3D x,std::vector <double> c)
{

    m_Pos = x;
    m_C = c;
    m_ID = id;
}


point::~point()
{
    
}


void point::UpdateBox(Vec3D *x)
{
m_pBox=x;
}
void point::UpdatePos(Vec3D x)
{
    m_Pos = x;
    
}
void point::UpdatepointUnitCell(UnitCell * z)
{
    m_PointUnitCell = z;
}
