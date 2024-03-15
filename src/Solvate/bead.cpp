#include "bead.h"

bead::bead(int id, const std::string& name, const std::string& type, const std::string& resname, int resid, double x, double y, double z)
    : m_X(x), m_Y(y), m_Z(z), m_ID(id), m_BeadName(name), m_BeadType(type), m_ResName(resname), m_Resid(resid), m_hasMol(false) {}

bead::bead(int id, const std::string& name, const std::string& type, const std::string& resname, int resid)
    : m_X(0.0), m_Y(0.0), m_Z(0.0), m_ID(id), m_BeadName(name), m_BeadType(type), m_ResName(resname), m_Resid(resid), m_hasMol(false) {}

bead::~bead() {}

void bead::UpdateBeadName(const std::string& s) {
    m_BeadName = s;
}

void bead::UpdateResName(const std::string& s) {
    m_ResName = s;
}

void bead::UpdateBox(Vec3D* x) {
    m_pBox = x;
}

void bead::UpdateXPos(double x) {
    m_X = x;
}

void bead::UpdateYPos(double y) {
    m_Y = y;
}

void bead::UpdateZPos(double z) {
    m_Z = z;
}

void bead::UpdatePos(Vec3D* B, double x, double y, double z) {
    m_pBox = B;
    UpdateXPos(x);
    UpdateYPos(y);
    UpdateZPos(z);
}

void bead::UpdatePos(double x, double y, double z) {
    UpdateXPos(x);
    UpdateYPos(y);
    UpdateZPos(z);
}

void bead::UpdateBeadUnitCell(UnitCell* z) {
    m_BeadUnitCell = z;
}

