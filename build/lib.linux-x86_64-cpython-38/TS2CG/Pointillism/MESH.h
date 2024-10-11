#if !defined(AFX_MESH_H_INCLUDED_)
#define AFX_MESH_H_INCLUDED_
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "CreateMashBluePrint.h"
#include "exclusion.h"


class MESH
{
public:

    MESH();
    ~MESH();

public:
    std::vector<vertex>         m_Vertex;
    std::vector<triangle>       m_Triangle;
    std::vector<links>          m_Links;
    std::vector<inclusion>      m_Inclusion;
    std::vector<exclusion>      m_Exclusion;
    Vec3D                       m_Box;

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
    std::vector<exclusion*>     m_pExclusion;

    void GenerateMesh(MeshBluePrint meshblueprint);
    MeshBluePrint Convert_Mesh_2_BluePrint(MESH *mesh);

};



#endif
