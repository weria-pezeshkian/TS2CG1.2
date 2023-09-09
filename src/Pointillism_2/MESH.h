#if !defined(AFX_MESH_H_INCLUDED_)
#define AFX_MESH_H_INCLUDED_
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "CreateMashBluePrint.h"
struct STRUC_Membrane_Parameters {
    double lambda;
    double kappa_geo;
    double kappa_normal;
    double kappa; // this is read by another parameters for now; later may change
    double kappa_g; // this is read by another parameters for now; later may change
};
class MESH
{
public:

    MESH();
    
    ~MESH();

    
private:
    std::vector<vertex>         m_Vertex;
    std::vector<triangle>       m_Triangle;
    std::vector<links>          m_Links;
    std::vector<inclusion>      m_Inclusion;
    Vec3D                       m_Box;
public:
    std::vector <InclusionType> m_InclusionType;
    std::vector <InclusionType*> m_pInclusionType;
    

    
    std::vector<vertex*>        m_pActiveV; // all the active vertices edge + surf
    std::vector<vertex*>        m_pSurfV; // all the active vertices  surf
    std::vector<vertex*>        m_pEdgeV;  // edge
    
    
    std::vector<links*>         m_pActiveL;   // all the links
    std::vector<links*>         m_pHL;
    std::vector<links*>         m_pMHL;
    std::vector<links*>         m_pEdgeL;
    
    
    std::vector<triangle*>      m_pActiveT;

    std::vector<inclusion*>     m_pInclusion;
    Vec3D                       *m_pBox;
    
    void GenerateMesh(MeshBluePrint meshblueprint, double kappa, double kappag, STRUC_Membrane_Parameters smp);
    MeshBluePrint Convert_Mesh_2_BluePrint(MESH *mesh);
    
    // MoveTrinagleFromAGroup2Another
    // MoveLinkFromAGroup2Another

    
    
};



#endif
