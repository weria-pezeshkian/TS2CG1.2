#if !defined(AFX_Solvate_H_SF4A21B8_C13C_1223_BF23_124095086234__INCLUDED_)
#define AFX_Solvate_H_SF4A21B8_C13C_1223_BF23_124095086234__INCLUDED_

#include "Def.h"
#include "bead.h"

class Solvate
{
public:
    
      Solvate(Argument *pArg);
	 ~Solvate();

private:
    void Bring2Box(std::vector<bead*> &Sysbead, Vec3D *Box);
    std::vector<bead> AddIons(std::vector<bead>& FullWaterBead, int Nposion, int Nnegion, const std::string& pName, const std::string& nName, int seed);

};

#endif
