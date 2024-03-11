#ifndef AFX_UnitCell_H_8P4B21B8_C13C_5648_BF23_444095086234__INCLUDED_
#define AFX_UnitCell_H_8P4B21B8_C13C_5648_BF23_444095086234__INCLUDED_

#include "Def.h"
#include "Vec3D.h"
#include <vector>

class bead;

class UnitCell {
public:
    UnitCell(int id, int i, int j, int k);
    ~UnitCell();

    inline const int GetID() const { return m_ID; }
    inline std::vector<bead*> GetBeadList() { return m_BeadList; }

public:
    void AddtoVertexList(bead* z);

private:
    int m_ID;
    std::vector<bead*> m_BeadList;
};

#endif
