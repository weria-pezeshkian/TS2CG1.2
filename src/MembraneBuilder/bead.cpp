

#include <stdio.h>
#include "bead.h"
//#include "links.h"
bead::bead(int id, std::string name, std::string type, std::string resname, int resid, double x, double y, double z)
{

m_X = x;
m_Y = y;
m_Z = z;
m_ID = id;
    m_BeadName = name;
    m_BeadType = type;
    m_ResName = resname;
    m_Resid = resid;

    m_hasMol =false;
}
bead::bead(int id, std::string name, std::string type, std::string resname, int resid)
{

m_ID=id;
m_BeadName = name;
m_BeadType = type;
m_X=0;
m_Y=0;
m_Z=0;
    m_ResName = resname;
    m_Resid = resid;
    m_hasMol =false;

}

bead::~bead() {
    
}


void bead::UpdateBox(Vec3D *x) {
m_pBox=x;
}
void bead::UpdateXPos(double x) {
		m_X=x;
}
void bead::UpdateYPos(double x){

		m_Y=x;
}
void bead::BringBeadInBox(Vec3D *pBox) {
    double boxX = (*pBox)(0);
    double boxY = (*pBox)(1);
    double boxZ = (*pBox)(2);

    // Adjust X position if it's outside the box bounds
    if (m_X >= boxX) {
        m_X = m_X - boxX;
    } else if (m_X < 0) {
        m_X = m_X + boxX;
    }
    
    if (m_Y >= boxY) {
        m_Y = m_Y - boxY;
    } else if (m_Y < 0) {
        m_Y = m_Y + boxY;
    }
    
    if (m_Z >= boxZ) {
        m_Z = m_Z - boxZ;
    } else if (m_Z < 0) {
        m_Z = m_Z + boxZ;
    }
    
}
void bead::UpdateZPos(double x)
{
    
    m_Z=x;

}
void bead::UpdatePos(Vec3D* B, double x, double y, double z)
{
    m_pBox=B;
    this->UpdateXPos(x);
    this->UpdateYPos(y);
    this->UpdateZPos(z);


}
void bead::UpdatePos(double x, double y, double z)
{
    this->UpdateXPos(x);
    this->UpdateYPos(y);
    this->UpdateZPos(z);
    
    
}
void bead::UpdateBeadUnitCell(UnitCell * z) {
    m_BeadUnitCell = z;
}
void bead::UpdateHasMol(bool z) {
    m_hasMol = z;
}












