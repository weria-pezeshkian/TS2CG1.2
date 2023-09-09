
#include <time.h>
#include <iomanip>
#include "Traj_XXX.h"


Traj_XXX::Traj_XXX(Vec3D *pBox)
{
 m_pBox=pBox;
m_Condition=true;
    Nfunction f;
    m_tsiPrecision = f.Int_to_String(18)+"."+f.Int_to_String(10);
}

Traj_XXX::~Traj_XXX()
{
    
}
/*void Traj_XXX::WriteTSI(int step ,  std::string filename , std::vector< vertex* > pver, std::vector< triangle* > ptriangle,  std::vector< inclusion* > pinc)
{
    FILE * output;
    output = fopen((filename).c_str(), "w");
    fprintf(output,"%5d\n",step);
    std::string format = "%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    fprintf(output,format.c_str(),(*m_pBox)(0),(*m_pBox)(1),(*m_pBox)(2));
    
    
    const char* ver="vertex";
    int size=pver.size();
    fprintf(output,"%s%20d\n",ver,size);
     format = "%5d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<vertex *>::iterator it = pver.begin() ; it != pver.end(); ++it)
    fprintf(output,format.c_str(),(*it)->GetVID(),(*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
    
    const char* tri="triangle";
    size = ptriangle.size();
    fprintf(output,"%s%20d\n",tri,size);
    for (std::vector<triangle *>::iterator it = ptriangle.begin() ; it != ptriangle.end(); ++it)
    fprintf(output,"%5d%5d%5d%5d\n",(*it)->GetTriID(),((*it)->GetV1())->GetVID(),((*it)->GetV2())->GetVID(),((*it)->GetV3())->GetVID());
    

    const char* inc="inclusion";
    size = pinc.size();
    fprintf(output,"%s%20d\n",inc,size);
    format = "%5d%5d%5d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<inclusion *>::iterator it = pinc.begin() ; it != pinc.end(); ++it)
    fprintf(output,format.c_str(),(*it)->GetID(),(*it)->GetTypeID(),((*it)->Getvertex())->GetVID(),((*it)->GetLDirection())(0),((*it)->GetLDirection())(1));
    
    fclose(output);

}*/
void Traj_XXX::WriteTSI(int step ,  std::string filename , std::vector< vertex* > pver, std::vector< triangle* > ptriangle,  std::vector< inclusion* > pinc , std::vector< exclusion* > pexc)
{
    FILE * output;
    output = fopen((filename).c_str(), "w");
    std::string format = "%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    const char* version="version 1.1";
    fprintf(output,"%s\n",version);
//------
    const char* box="box";
    fprintf(output,"%s%18.10lf%18.10lf%18.10lf\n",box,(*m_pBox)(0),(*m_pBox)(1),(*m_pBox)(2));

    const char* ver="vertex";
    int size=pver.size();
    fprintf(output,"%s%20d\n",ver,size);
    format = "%5d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<vertex *>::iterator it = pver.begin() ; it != pver.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetVID(),(*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
    
    const char* tri="triangle";
    size = ptriangle.size();
    fprintf(output,"%s%20d\n",tri,size);
    for (std::vector<triangle *>::iterator it = ptriangle.begin() ; it != ptriangle.end(); ++it)
        fprintf(output,"%10d%10d%10d%10d\n",(*it)->GetTriID(),((*it)->GetV1())->GetVID(),((*it)->GetV2())->GetVID(),((*it)->GetV3())->GetVID());
    
    
    const char* inc="inclusion";
    size = pinc.size();
    fprintf(output,"%s%20d\n",inc,size);
    format = "%10d%10d%10d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<inclusion *>::iterator it = pinc.begin() ; it != pinc.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetID(),(*it)->GetTypeID(),((*it)->Getvertex())->GetVID(),((*it)->GetLDirection())(0),((*it)->GetLDirection())(1));
    
    if(pexc.size()!=0)
    {
        const char* exc="exclusion";
        size = pexc.size();
        fprintf(output,"%s%20d\n",exc,size);
        format = "%10d%10d%"+m_tsiPrecision+"lf\n";
        for (std::vector<exclusion *>::iterator it = pexc.begin() ; it != pexc.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetID(),((*it)->Getvertex())->GetVID(),(*it)->GetRadius());
    }
    
    
    
    
    fclose(output);
    
}
void Traj_XXX::ReadTSI(std::string filename)
{
    
    std::ifstream input;
    input.open(filename.c_str());
    std::string version = "0.0";
    std::string str;
    Nfunction f;


    while (true)
    {
        getline(input,str);
        if(input.eof())
            break;
        std::vector<std::string> S = f.split(str);
    
        if(S.size()>1)
        if(S.at(0)=="version")
        {
            version = S.at(1);
        }
    }
    input.close();
    if(version=="0.0")
    {
        ReadTSI1(filename);
    }
    else
    {
        ReadTSI2(filename);
        std::cout<<"---> Note: the tsi file version is "<<version<<"\n";

    }
}
    
void Traj_XXX::ReadTSI2(std::string filename)
{
    
    //=======================
    //===== Clear the storage
    m_Vertex.clear();
    m_Triangle.clear();
    m_Links.clear();
    m_Inclusion.clear();
    m_pAllV.clear();
    m_pAllT.clear();
    m_pLinks.clear();
    m_pHL.clear();
    m_pMHL.clear();
    m_pInclusion.clear();
    //======================
    int step,nver,ntr,ninc,id,tid,vid;
    double x,y,z;
    int v1,v2,v3,domain;
    std::string version,str;
    Nfunction f;
    std::ifstream input;
    input.open(filename.c_str());
    (*m_pBox)(0) = 100;
    (*m_pBox)(1) = 100;
    (*m_pBox)(2) = 100;
    while (true)
    {
        input>>str;
        if(input.eof())
            break;
        
        if(str=="box")
        {
            getline(input,str);
            std::vector<std::string> S = f.split(str);
            if(S.size()<3)
            {
                std::cout<<"error ---> information of the box is not sufficent in the tsi file \n";
                break;
            }
            else
            {
                (*m_pBox)(0) = f.String_to_Double(S.at(0));
                (*m_pBox)(1) = f.String_to_Double(S.at(1));
                (*m_pBox)(2) = f.String_to_Double(S.at(2));
            }
        }
        else if(str=="vertex")
        {
            input>>nver;
            double x,y,z;
            int id,domain;
            getline(input,str);

            for (int i=0;i<nver;i++)
            {
                getline(input,str);
                std::vector<std::string> S = f.split(str);

                if(S.size()<4)
                {
                    std::cout<<"error ---> information of the vertex "<<i<<" is not sufficent in the tsi file \n";
                    m_Condition =false;
                    std::cout<<str<<" \n";

                    break;
                }
                else
                {
                    id=f.String_to_Int(S.at(0));
                    x=f.String_to_Double(S.at(1));
                    y=f.String_to_Double(S.at(2));
                    z=f.String_to_Double(S.at(3));
                    vertex v(id,x,y,z);
                    v.UpdateBox(m_pBox);
                    v.UpdateClusterID(1);
                    if(S.size()>4)
                    {
                        domain=f.String_to_Double(S.at(4));
                        v.UpdateDomainID(domain);
                    }
                    m_Vertex.push_back(v);
                }
            }

        }
        else if(str=="triangle")
        {
            input>>ntr;
            getline(input,str);

            int id,v1,v2,v3;

            for (int i=0;i<ntr;i++)
            {
                getline(input,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<4)
                {
                    std::cout<<"error ---> information of the trinagles  "<<i<<" is not sufficent in the tsi file \n";
                    m_Condition =false;

                    std::cout<<str<<" \n";

                    break;
                }
                else
                {
                    id=f.String_to_Int(S.at(0));
                    v1=f.String_to_Int(S.at(1));
                    v2=f.String_to_Int(S.at(2));
                    v3=f.String_to_Int(S.at(3));

                    triangle T(id,&(m_Vertex.at(v1)),&(m_Vertex.at(v2)),&(m_Vertex.at(v3)));
                    m_Triangle.push_back(T);
                    
                }
            }
        }
        else if(str=="inclusion")
        {
            input>>ninc;
            getline(input,str);

            int id,vid,type;
            double Dx, Dy;
            for (int i=0;i<ninc;i++)
            {
                getline(input,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<5)
                {
                    std::cout<<"error ---> information of the inclusion "<<i<<" is not sufficent in the tsi file \n";
                    m_Condition =false;

                    std::cout<<str<<" \n";

                    break;
                }
                else
                {
                id=f.String_to_Int(S.at(0));
                type=f.String_to_Int(S.at(1));
                vid=f.String_to_Int(S.at(2));
                Dx=f.String_to_Double(S.at(3));
                Dy=f.String_to_Double(S.at(4));
                double norm = sqrt(Dx*Dx+Dy*Dy);
                    Dy=Dy/norm;
                    Dx=Dx/norm;
                    
                inclusion Tinc(id);
                Tinc.Updatevertex(&(m_Vertex.at(vid)));
                Vec3D D;
                D(0)=Dx;D(1)=Dy;D(2)=0;
                Tinc.UpdateType("P"+f.Int_to_String(type),type);
                Tinc.UpdateLocalDirection(D);
                m_Inclusion.push_back(Tinc);
                (m_Vertex.at(vid)).UpdateOwnInclusion(true);
                }

            }

        }
        else if(str=="exclusion")
        {
            input>>ninc;
            getline(input,str);
            
            int id,vid;
            double r;
            for (int i=0;i<ninc;i++)
            {
                getline(input,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<3)
                {
                    std::cout<<"error ---> information of the inclusion "<<i<<" is not sufficent in the tsi file \n";
                    m_Condition =false;
                    
                    std::cout<<str<<" \n";
                    
                    break;
                }
                else
                {
                    id=f.String_to_Int(S.at(0));
                    vid=f.String_to_Int(S.at(1));
                    r=f.String_to_Double(S.at(2));
          
                    
                    exclusion Tinc(id);
                    Tinc.Updatevertex(&(m_Vertex.at(vid)));
                    Tinc.UpdateRadius(r);
                    m_Exclusion.push_back(Tinc);
                }
                
            }
            
        }
        else
        {
            getline(input,str);

        }
        

    }
    input.close();
    
    
    for (std::vector<inclusion>::iterator it = m_Inclusion.begin() ; it != m_Inclusion.end(); ++it)
    m_pInclusion.push_back(&(*it));
    
    
    for (std::vector<exclusion>::iterator it = m_Exclusion.begin() ; it != m_Exclusion.end(); ++it)
        m_pExclusion.push_back(&(*it));

    //======= we finished the reading.
    for (std::vector<vertex>::iterator it = m_Vertex.begin() ; it != m_Vertex.end(); ++it)
        m_pAllV.push_back(&(*it));

    int li=-1;
    
    for (std::vector<triangle>::iterator it = m_Triangle.begin() ; it != m_Triangle.end(); ++it)
    {
        
        m_pAllT.push_back(&(*it));
        (it->GetV1())->AddtoTraingleList(&(*it));
        (it->GetV1())->AddtoNeighbourVertex((it->GetV2()));
        (it->GetV2())->AddtoTraingleList(&(*it));
        (it->GetV2())->AddtoNeighbourVertex((it->GetV3()));
        (it->GetV3())->AddtoTraingleList(&(*it));
        (it->GetV3())->AddtoNeighbourVertex((it->GetV1()));
        
        /// create links
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        
        links l1(id1,it->GetV1(),it->GetV2(),&(*it));
        l1.UpdateV3(it->GetV3());
        
        links l2(id2,it->GetV2(),it->GetV3(),&(*it));
        l2.UpdateV3(it->GetV1());
        
        links l3(id3,it->GetV3(),it->GetV1(),&(*it));
        l3.UpdateV3(it->GetV2());
        m_Links.push_back(l1);
        m_Links.push_back(l2);
        m_Links.push_back(l3);
        
    }
     li=-1;
    for (std::vector<triangle>::iterator it = m_Triangle.begin() ; it != m_Triangle.end(); ++it)
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
        
        
        (it->GetV1())->AddtoLinkList(l1);
        (it->GetV2())->AddtoLinkList(l2);
        (it->GetV3())->AddtoLinkList(l3);

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
            std::string A=" Warning : This system does not look like a closed system or the triangles orination are not consistent ";
            std::cout<<A<<"\n";
            f.Write_One_LogMessage(A);
        }
        
    }
    int nomirror=0;
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        if(it->GetMirrorFlag()==false)
        {
            nomirror++;
        }
        
    }
    if(nomirror!=0)
    {
        std::string A=" Warning : there is links without mirror ";
        std::cout<<nomirror<<" links without mirror \n";
        f.Write_One_LogMessage(A);
    }
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        
        m_pLinks.push_back(&(*it));
    }
    
    for (std::vector<inclusion*>::iterator it = m_pInclusion.begin() ; it != m_pInclusion.end(); ++it)
    {
        ((*it)->Getvertex())->UpdateInclusion(*it);
    }


}
void Traj_XXX::ReadTSI1(std::string filename)
{
    //=======================
    //===== Clear the storage
    m_Vertex.clear();
    m_Triangle.clear();
    m_Links.clear();
    m_Inclusion.clear();
    m_pAllV.clear();
    m_pAllT.clear();
    m_pLinks.clear();
    m_pHL.clear();
    m_pMHL.clear();
    m_pInclusion.clear();
    //======================
    Nfunction f;
    char str[10];
    int step,nver,ntr,ninc,id,tid,vid;
    double x,y,z;
    
    int v1,v2,v3;
    FILE * input;
    input = fopen(filename.c_str(), "r");
    int rf = fscanf(input, "%d" ,&step );
    rf = fscanf(input, "%lf%lf%lf" ,&x,&y,&z );
    (*m_pBox)(0) = x;
    (*m_pBox)(1) = y;
    (*m_pBox)(2) = z;
    
    rf = fscanf(input, "%s%d" ,str,&nver );
    
    
    
    for (int i=0;i<nver;i++)
    {
        rf = fscanf(input, "%d%lf%lf%lf" ,&id,&x,&y,&z);
        vertex v(id,x,y,z);
        v.UpdateBox(m_pBox);
        v.UpdateClusterID(1);
        m_Vertex.push_back(v);
    }
    rf = fscanf(input, "%s%d" ,str,&ntr );
    
    for (int i=0;i<ntr;i++)
    {
        rf = fscanf(input, "%d%d%d%d" ,&id,&v1,&v2,&v3);
        bool pr=true;
        triangle T(id,&(m_Vertex.at(v1)),&(m_Vertex.at(v2)),&(m_Vertex.at(v3)));
        
        m_Triangle.push_back(T);
    }
    rf = fscanf(input, "%s%d" ,str,&ninc );
    for (int i=0;i<ninc;i++)
    {
        Vec3D D;
        rf = fscanf(input, "%d%d%d%lf%lf" ,&id,&tid,&vid,&x,&y);
        D(0)=x;D(1)=y;D(2)=0;
        inclusion Tinc(id);
        Tinc.Updatevertex(&(m_Vertex.at(vid)));
        Tinc.UpdateType("P"+f.Int_to_String(tid),tid);
        Tinc.UpdateLocalDirection(D);
        m_Inclusion.push_back(Tinc);
        (m_Vertex.at(vid)).UpdateOwnInclusion(true);
        
        
        
    }
    fclose(input);
    for (std::vector<inclusion>::iterator it = m_Inclusion.begin() ; it != m_Inclusion.end(); ++it)
        m_pInclusion.push_back(&(*it));
    //======= we finished the reading.
    for (std::vector<vertex>::iterator it = m_Vertex.begin() ; it != m_Vertex.end(); ++it)
        m_pAllV.push_back(&(*it));
    
    int li=-1;
    
    for (std::vector<triangle>::iterator it = m_Triangle.begin() ; it != m_Triangle.end(); ++it)
    {
        
        m_pAllT.push_back(&(*it));
        (it->GetV1())->AddtoTraingleList(&(*it));
        (it->GetV1())->AddtoNeighbourVertex((it->GetV2()));
        (it->GetV2())->AddtoTraingleList(&(*it));
        (it->GetV2())->AddtoNeighbourVertex((it->GetV3()));
        (it->GetV3())->AddtoTraingleList(&(*it));
        (it->GetV3())->AddtoNeighbourVertex((it->GetV1()));
        
        /// create links
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        
        links l1(id1,it->GetV1(),it->GetV2(),&(*it));
        l1.UpdateV3(it->GetV3());
        
        links l2(id2,it->GetV2(),it->GetV3(),&(*it));
        l2.UpdateV3(it->GetV1());
        
        links l3(id3,it->GetV3(),it->GetV1(),&(*it));
        l3.UpdateV3(it->GetV2());
        m_Links.push_back(l1);
        m_Links.push_back(l2);
        m_Links.push_back(l3);
        
    }
    li=-1;
    for (std::vector<triangle>::iterator it = m_Triangle.begin() ; it != m_Triangle.end(); ++it)
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
        
        
        (it->GetV1())->AddtoLinkList(l1);
        (it->GetV2())->AddtoLinkList(l2);
        (it->GetV3())->AddtoLinkList(l3);
        
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
            std::string A=" Warning : This system does not look like a closed system or the triangles orination are not consistent ";
            std::cout<<A<<"\n";
            f.Write_One_LogMessage(A);
        }
        
    }
    int nomirror=0;
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        if(it->GetMirrorFlag()==false)
        {
            nomirror++;
        }
        
    }
    if(nomirror!=0)
    {
        std::string A=" Warning : there is links without mirror ";
        std::cout<<nomirror<<" links without mirror \n";
        f.Write_One_LogMessage(A);
    }
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        
        m_pLinks.push_back(&(*it));
    }
    
    for (std::vector<inclusion*>::iterator it = m_pInclusion.begin() ; it != m_pInclusion.end(); ++it)
    {
        ((*it)->Getvertex())->UpdateInclusion(*it);
    }
    
}


