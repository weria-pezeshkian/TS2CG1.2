#if !defined(AFX_BackMap_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_BackMap_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

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
#include "ReadDTSFolder.h"
#include "Data_Structure.h"



class BackMap
{
public:
    
	BackMap(Argument *pArgu);
	virtual ~BackMap();
    


public:

private:
    std::vector<bead*> m_pAllBeads;
    Vec3D *m_pBox;
    Vec3D  m_Box;
    std::vector<bead> m_FinalBeads;
    std::string m_FinalOutputGroFileName;
    std::string m_FinalTopologyFileName;

    std::map<std::string , MolType>  m_MoleculesType;
    int m_ResID;
    std::vector<inclusion>  m_RandomInc;
    std::vector<inclusion*>  m_pInc;
    std::vector<ExcludedVolumeBeads>  m_ExcludeBeads;
    double m_Iter;
    std::string m_InclusionDirectionType;
    std::vector<exclusion*>  m_pExc;

private:
    //=== since 2024
    void ExcludePointsUsingExclusion(std::vector<exclusion*>&, std::vector<point*>&, std::vector<point*>&);
    bool CheckProteinInfo (std::map<int , ProteinList>&, std::map<std::string , MolType>&, std::vector<inclusion*> &);
    
    
    
    //======== old functions
    double dist2between2Points(Vec3D X1,Vec3D X2);
    void WriteFinalGroFile();
    bool AnyBeadWithinCutoff(UnitCell *,Vec3D Pos);
    Tensor2 Rz(double cos, double sin);
    Tensor2 TransferMatLG(Vec3D N, Vec3D t1, Vec3D t2);
    void CreateRandomInclusion();
    std::map<int , ProteinList>  m_TotalProteinList;
    bool m_Renormalizedlipidratio;
    double    m_TotalAreaUp ;
    double    m_TotalAreaDown ;
    double    m_AvailPointUp ;
    double    m_AvailPointDown ;
    bool m_monolayer;
    std::vector<point*>  m_point1;
    std::vector<point*>  m_point2;

   // void CreateWallBead(std::vector<point*>  p1, std::vector<point*>  p2);
    bool FindProteinList(std::string str);
    void GenProtein(MolType moltype, int , Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2);
    void GenLipid(MolType moltype, int , Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2);

};


#endif
