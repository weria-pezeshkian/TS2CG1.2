#include "MESH.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
MESH class for quick access 
 */
MESH::MESH()
{
    
    
}
MESH::~MESH()
{
    
}
void MESH::GenerateMesh(MeshBluePrint meshblueprint, double kappa, double kappag, STRUC_Membrane_Parameters smp)
{
    m_Box = meshblueprint.simbox;
    m_pBox = &m_Box;
    m_InclusionType = meshblueprint.binctype;
    
    for (std::vector<InclusionType>::iterator it = m_InclusionType.begin() ; it != m_InclusionType.end(); ++it)
        m_pInclusionType.push_back(&(*it));
    
    
    // Making vertices
    for (std::vector<Vertex_Map>::iterator it = (meshblueprint.bvertex).begin() ; it != (meshblueprint.bvertex).end(); ++it)
    {
            vertex v(it->id,it->x,it->y,it->z);
            v.UpdateBox(m_pBox);
            v.UpdateGroup(it->domain);
            v.UpdateKappa(kappa/2.0,kappag);
            v.m_Lambda = smp.lambda;
            v.m_KGC = smp.kappa_geo;
            v.m_KNC = smp.kappa_normal;
            m_Vertex.push_back(v);
    }
//===== Make exclution [since June, 2023]
       
       std::vector<int> excluded_ver = meshblueprint.excluded_id; // this vertices should be excluded
       for (std::vector<int>::iterator it = excluded_ver.begin() ; it != excluded_ver.end(); ++it)
       {
           ((meshblueprint.bvertex).at((*it))).include = false;
       }
       int t=0;
       for (std::vector<vertex>::iterator it = m_Vertex.begin() ; it != m_Vertex.end(); ++it)
       {
               if(((meshblueprint.bvertex).at((t))).include == true)
                   m_pActiveV.push_back(&(*it));
           t++;
       }
    //== since the vertex number remain constant, it is better to change the id of active vertices
    t=0;
    for (std::vector<vertex*>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it)
    {
        (*it)->UpdateVID(t);
        t++;
    }
    // Making triangles
    t=0;
    int temid = 0;
    for (std::vector<Triangle_Map>::iterator it = (meshblueprint.btriangle).begin() ; it != (meshblueprint.btriangle).end(); ++it)
    {
        bool pr=true;
        
        int vid1 = ((meshblueprint.btriangle).at(t)).v1;
        int vid2 = ((meshblueprint.btriangle).at(t)).v2;
        int vid3 = ((meshblueprint.btriangle).at(t)).v3;

        if(((meshblueprint.bvertex).at(vid1)).include == true)
        if(((meshblueprint.bvertex).at(vid2)).include == true)
        if(((meshblueprint.bvertex).at(vid3)).include == true)
        {
        triangle T(temid,&(m_Vertex.at(it->v1)),&(m_Vertex.at(it->v2)),&(m_Vertex.at(it->v3)));
        m_Triangle.push_back(T);
            temid++;
        }
        t++;
    }
    //make inclusions
    for (std::vector<Inclusion_Map>::iterator it = (meshblueprint.binclusion).begin() ; it != (meshblueprint.binclusion).end(); ++it)
    {
        inclusion Tinc(it->id);
        if(m_Vertex.size()<it->vid+1)
        {
        std::cout<<"----> Error: Inclusion vertex id is out of range "<<std::endl;
            exit(0);
        }
        Tinc.Updatevertex(&(m_Vertex.at(it->vid)));
        Tinc.UpdateInclusionTypeID(it->tid);
        Vec3D D(it->x,it->y,0);
        Tinc.UpdateLocalDirection(D);
        m_Inclusion.push_back(Tinc);
        (m_Vertex.at(it->vid)).UpdateOwnInclusion(true);
    }
    for (std::vector<inclusion>::iterator it = m_Inclusion.begin() ; it != m_Inclusion.end(); ++it)
        m_pInclusion.push_back(&(*it));


    
    t=0;
    for (std::vector<triangle>::iterator it = m_Triangle.begin() ; it != m_Triangle.end(); ++it)
        m_pActiveT.push_back(&(*it));

// end Make exclution
    
    
    int li=-1;
    
    for (std::vector<triangle*>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it)
    {
        ((*it)->GetV1())->AddtoTraingleList((*it));
        ((*it)->GetV1())->AddtoNeighbourVertex(((*it)->GetV2()));
        ((*it)->GetV2())->AddtoTraingleList((*it));
        ((*it)->GetV2())->AddtoNeighbourVertex(((*it)->GetV3()));
        ((*it)->GetV3())->AddtoTraingleList((*it));
        ((*it)->GetV3())->AddtoNeighbourVertex(((*it)->GetV1()));
        
        /// create links
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        
        links l1(id1,(*it)->GetV1(),(*it)->GetV2(),(*it));
        l1.UpdateV3((*it)->GetV3());
        
        links l2(id2,(*it)->GetV2(),(*it)->GetV3(),(*it));
        l2.UpdateV3((*it)->GetV1());
        
        links l3(id3,(*it)->GetV3(),(*it)->GetV1(),(*it));
        l3.UpdateV3((*it)->GetV2());
        m_Links.push_back(l1);
        m_Links.push_back(l2);
        m_Links.push_back(l3);
        
    }
    li=-1;
    for (std::vector<triangle*>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it)
    {
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        links * l1=&(m_Links.at(id1));
        links * l2=&(m_Links.at(id2));
        links * l3=&(m_Links.at(id3));
        l1->UpdateNeighborLink1(l2);
        l1->UpdateNeighborLink2(l3);
        l2->UpdateNeighborLink1(l3);
        l2->UpdateNeighborLink2(l1);
        l3->UpdateNeighborLink1(l1);
        l3->UpdateNeighborLink2(l2);
        
        
        ((*it)->GetV1())->AddtoLinkList(l1);
        ((*it)->GetV2())->AddtoLinkList(l2);
        ((*it)->GetV3())->AddtoLinkList(l3);
        
    }
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        bool foundM=false;
        if((it)->GetMirrorFlag()==true)
        {
            m_pMHL.push_back(it->GetMirrorLink());
            m_pHL.push_back(&(*it));
            foundM = true;
        }
        else
        {
            vertex *v1=it->GetV1();
            vertex *v2=it->GetV2();
            
            std::vector<links*>  lList = v2->GetVLinkList();
            for (std::vector<links*>::iterator it2 = lList.begin() ; it2 != lList.end(); ++it2)
            {
                if(((*it2)->GetV2())->GetVID()==v1->GetVID())
                {
                    it->UpdateMirrorLink((*it2));
                    (*it2)->UpdateMirrorLink(&(*it));
                    it->UpdateMirrorFlag(true);
                    (*it2)->UpdateMirrorFlag(true);
                    foundM = true;
                    break;
                }
            }
        }
        if(foundM == false)
        {
            m_pEdgeL.push_back(&(*it));
            it->m_LinkType = 1;
        }
        
    }
    int edgelink=0;
    for (std::vector<links*>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it)
        edgelink++;

    if(edgelink!=0)
    {
        std::cout<<"----> Note: "<<edgelink<<" links at the edge \n";
        std::cout<<"----> Note (Warning): the system is not closed! \n";

    }
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        
        m_pActiveL.push_back(&(*it));
    }
    //==== Getting the edge vertex from link
    for (std::vector<links*>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it)
    {
        m_pEdgeV.push_back((*it)->GetV1());
        ((*it)->GetV1())->m_pEdgeLink = *it;
        ((*it)->GetV2())->m_pPrecedingEdgeLink = *it;
        ((*it)->GetV2())->AddtoNeighbourVertex((*it)->GetV1());
        ((*it)->GetV1())->m_VertexType = 1;
    }
    for (std::vector<vertex*>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it)
    {
        bool isedge = false;
        for (std::vector<vertex*>::iterator it2 = m_pEdgeV.begin() ; it2 != m_pEdgeV.end(); ++it2)
        {
            if((*it2)->GetVID()==(*it)->GetVID())
                isedge = true;
        }
        if(isedge==false)
        m_pSurfV.push_back(*it);
    }
    
    //
    
    
    for (std::vector<inclusion*>::iterator it = m_pInclusion.begin() ; it != m_pInclusion.end(); ++it)
    {
        ((*it)->Getvertex())->UpdateInclusion((*it));
        int inc_typeid=(*it)->GetInclusionTypeID();
        if(m_InclusionType.size()-1<inc_typeid)
        {
            std::cout<<" Error: inclusion with typeid of "<<inc_typeid<<" has not been defined \n";
            exit(0);
        }
        (*it)->UpdateInclusionType(&(m_InclusionType.at(inc_typeid)));
    }
    // =======
    // ==== info of the mesh
    
    std::cout<<"---> active vertex "<<m_pActiveV.size()<<" surf vertex "<<m_pSurfV.size()<<"  edge vertex "<<m_pEdgeV.size()<<" -- \n";




    
    //WritevtuFiles VTU(pState);
    //std::string file="ini_Mesh.vtu";
    //VTU.Writevtu(m_pAllV,m_pAllT,m_pAllLinks,file);
}
//===========================================================
// Note, the converted blue print will not have the exclusions
//
MeshBluePrint MESH::Convert_Mesh_2_BluePrint(MESH *mesh)
{
    MeshBluePrint BluePrint;
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    std::vector <InclusionType> binctype;  // a vector containing all inclsuion type and a default one
    Vec3D simbox;
    
    // vertex member of the blue print
    std::vector<vertex*> pV = mesh->m_pActiveV;
    for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        Vertex_Map tvm;
        tvm.x = (*it)->GetVXPos();
        tvm.y = (*it)->GetVYPos();
        tvm.z = (*it)->GetVZPos();
        tvm.id = (*it)->GetVID();
        tvm.domain = (*it)->GetGroup();
        bvertex.push_back(tvm);
    }
    // triangle map member of the blue print
    std::vector<triangle*> pT = mesh->m_pActiveT;
    for (std::vector<triangle *>::iterator it = pT.begin() ; it != pT.end(); ++it)
    {
        Triangle_Map ttm;
        ttm.v1 = ((*it)->GetV1())->GetVID();
        ttm.v2 = ((*it)->GetV2())->GetVID();
        ttm.v3 = ((*it)->GetV3())->GetVID();
        ttm.id = (*it)->GetTriID();
        btriangle.push_back(ttm);
    }
    
    // inclusion map member of the blue print
    std::vector<inclusion*> pInc = mesh->m_pInclusion;
    for (std::vector<inclusion *>::iterator it = pInc.begin() ; it != pInc.end(); ++it)
    {
        Inclusion_Map tim;
        tim.x = ((*it)->GetLDirection())(0);
        tim.y = ((*it)->GetLDirection())(1);
        tim.vid = ((*it)->Getvertex())->GetVID();
        tim.tid = ((*it)->GetInclusionType())->ITid;
        tim.id = ((*it)->GetID());
        binclusion.push_back(tim);

    }
    // inclusion type map member of the blue print
    BluePrint.binctype = mesh->m_InclusionType;
    // Add other map into the mesh map
    BluePrint.bvertex = bvertex;
    BluePrint.btriangle = btriangle;
    BluePrint.binclusion = binclusion;
    BluePrint.simbox = *(mesh->m_pBox);
    
    return BluePrint;
}

