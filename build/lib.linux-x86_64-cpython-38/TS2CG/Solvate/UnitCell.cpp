#include <stdio.h>
#include "UnitCell.h"

UnitCell::UnitCell(int id, int i, int j, int k) {
    m_ID = id;
}

UnitCell::~UnitCell() {
    
}

void UnitCell::AddtoVertexList(bead* z) {
    m_BeadList.push_back(z);
}
