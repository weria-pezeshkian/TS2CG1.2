

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
    m_Pos(0) = x;
    m_Pos(1) = y;
    m_Pos(2) = z;

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
    
    m_Pos(0) = 0;
    m_Pos(1) = 0;
    m_Pos(2) = 0;
    m_C.push_back(0);
    m_C.push_back(0);

}

bead::~bead()
{
    
}

void  bead::UpdateC(std::vector <double> c)
{
    m_C = c;
}

void bead::UpdateBox(Vec3D *x)
{
m_pBox=x;
}
void bead::UpdateXPos(double x)
{
	if(x>=(*m_pBox)(0))
	{
		m_X=x-(*m_pBox)(0);
	}
	else if(x<0)
	{
		m_X=x+(*m_pBox)(0);
	}
	else
	{
		m_X=x;
	}
    
    m_Pos(0) = m_X;
    
}
void bead::UpdateYPos(double x)
{
	if(x>=(*m_pBox)(1))
	{
		m_Y=x-(*m_pBox)(1);
	}
	else if(x<0)
	{
		m_Y=x+(*m_pBox)(1);
	}
	else
	{
		m_Y=x;
	}
    
    m_Pos(1) = m_Y;

}
void bead::UpdateZPos(double x)
{
	if(x>=(*m_pBox)(2))
	{
		m_Z=x-(*m_pBox)(2);
	}
	else if(x<0)
	{
		m_Z=x+(*m_pBox)(2);
	}
	else
	{
		m_Z=x;
	}
    
    m_Pos(2) = m_Z;

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
void bead::UpdateBeadUnitCell(UnitCell * z)
{
    m_BeadUnitCell = z;
}
void bead::UpdateBeadMol(molecules * z)
{
    m_pBeadMol = z;
}
void bead::UpdateHasMol(bool z)
{
    m_hasMol = z;
}












