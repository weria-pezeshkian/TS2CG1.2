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
#include "GenDomains.h"



class BackMap
{
public:
    
	BackMap(Argument *pArgu);
	virtual ~BackMap();
    


public:

    
private:
    // final data
    std::string m_FinalOutputGroFileName;
    std::string m_FinalTopologyFileName;
    
private:
    // all the maps
    std::map<std::string , MolType>  m_map_MolName2MoleculesType;
    std::map<int , ProteinList>  m_map_IncID2ProteinLists;
    std::vector<bead> m_FinalBeads;         // all the beads generated at the end

private:

    std::vector<bead*> m_pAllBeads;
    Vec3D *m_pBox;
    Vec3D  m_Box;




    int m_ResID;
    std::vector<ExcludedVolumeBeads>  m_ExcludeBeads;
    double m_Iter;
    std::string m_InclusionDirectionType;


private:
    //=== since 2024
    void ExcludePointsUsingExclusion(std::vector<exclusion*>&, std::vector<point*>&, std::vector<point*>&);
    bool CheckProteinInfo (std::map<int , ProteinList>&, std::map<std::string , MolType>&, std::vector<inclusion*> &);
    bool PlaceProteins(std::vector<point*> &PointUp, std::vector<inclusion*>  &pInc);
    bool RemovePointsCloseToBeadList(std::vector<point*> &PointUp, std::vector<point*> &PointDown, std::vector<bead*> vpbeads, double RCutOff, Vec3D* m_pBox);
    bool GenLipidsForADomain(Domain *pdomain); // generates all the lipid for a specific domain
    bool GenTopologyFile(std::vector<Domain*>, int wbeadno); // generates topology file

    
    //======== old functions
    double dist2between2Points(Vec3D X1,Vec3D X2);
    void WriteFinalGroFile();
    bool AnyBeadWithinCutoff(UnitCell *,Vec3D Pos);
    Tensor2 Rz(double cos, double sin);
    Tensor2 TransferMatLG(Vec3D N, Vec3D t1, Vec3D t2);
    std::vector<inclusion> CreateRandomInclusion(std::vector<point*> &PointUp);
    //bool m_Renormalizedlipidratio;
    bool m_monolayer;

   // void CreateWallBead(std::vector<point*>  p1, std::vector<point*>  p2);
    bool FindProteinList(std::string str);
    void GenProtein(MolType moltype, int , Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2);
    void GenLipid(MolType moltype, int , Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2);
    void Welldone();

};


#endif
