

#include <stdio.h>
#include "ReadDTSFolder.h"
ReadDTSFolder::ReadDTSFolder()
{
    
}
void ReadDTSFolder::Read(std::string foldername)
{



    std::string file1 = "./"+foldername+"/OuterBM.dat";
    std::string file2 = "./"+foldername+"/InnerBM.dat";
    std::string file3 = "./"+foldername+"/IncData.dat";
    std::string file4 = "./"+foldername+"/ExcData.dat";

    bool health = true;
    bool monolayer =false;
    if(FileExist(file1) == false)
    {
        std::cout<<"--> error: file with file name "<<file1 <<" does not exist \n";
        health = false;
    }
    if(FileExist(file2) == false)
    {
        monolayer = true;
    }
    if(FileExist(file3) == false)
    {
        std::cout<<"--> no inclsuion file is provided, we will generate a random distribution of proteins if information is provided in STR file \n";
    }
    
    m_OuterPoint = ReadPointObjects(file1,1);
    if(monolayer==false)
    m_InnerPoint = ReadPointObjects(file2,-1);
    
    if(FileExist(file3) == true)
    {
        std::cout<<"--> inclusion file is provided, we will generate proteins according to this file \n";
        m_Inclusion =  ReadInclusionObjects(file3);

    }
    if(FileExist(file4) == true)
    {
        std::cout<<"--> exclusion file is provided, meaning the system contains pores \n";
        m_Exclusion =  ReadExclusionObjects(file4);
    }
    


}
ReadDTSFolder::~ReadDTSFolder()
{
    
}



std::vector<point> ReadDTSFolder::ReadPointObjects(std::string file, int lay)
{
    

   //  char str = new str[1000];
    char str1[256];
    char str2[256];

    FILE *fdtspoins;
    fdtspoins = fopen(file.c_str(), "r");
    char layer[256];
    float Lx,Ly,Lz;
    

    int NoPoints;

    if (fdtspoins == NULL){
        printf(" Error: Could not open file %s",file.c_str());
    }

    if(lay==1)
    int readafile =  fscanf(fdtspoins,"%s%f%f%f",str2,&Lx,&Ly,&Lz);

    int readafile = fscanf(fdtspoins,"%s%s%s%d%s",str2,str2,str2,&NoPoints,str2);

    m_Box(0) = Lx;
    m_Box(1) = Ly;
    m_Box(2) = Lz;
    
        
    
    bool check = fgets(str1, sizeof(str1), fdtspoins);
    check = fgets(str1, sizeof(str1), fdtspoins);
    check = fgets(layer, sizeof(layer), fdtspoins);
    std::string LayerName=layer;
    LayerName.pop_back();
    LayerName.pop_back();
    LayerName.pop_back();
    LayerName.erase(LayerName.begin());
    LayerName.erase(LayerName.begin());

    std::vector<point>  AllPoint;


    float area,x,y,z,nx,ny,nz,p1x,p1y,p1z,p2x,p2y,p2z,c1,c2;
    int id,domainid;
    for (int i=0;i<NoPoints;i++)
    {
        check = fscanf(fdtspoins,"%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&id,&domainid,&area,&x,&y,&z,&nx,&ny,&nz,&p1x,&p1y,&p1z,&p2x,&p2y,&p2z,&c1,&c2);
        

        if(area==0)
        {
            std::cout<<"point id "<<id<<"  has zero area \n";
        }

        Vec3D X(x,y,z);
        Vec3D N(nx,ny,nz);
        Vec3D P1(p1x,p1y,p1z);
        Vec3D P2(p2x,p2y,p2z);
        std::vector <double> C;
        C.push_back(c1);
        C.push_back(c2);
        point p(id, area,X,N, P1, P2, C );
        if(lay==-1)
        p.UpdateUpperLayer(false);
        p.UpdateDomainID(domainid);
        AllPoint.push_back(p);
    }
    
    return AllPoint;

}
std::vector<inclusion> ReadDTSFolder::ReadInclusionObjects(std::string file)
{
//    inclusion(int id, int typeID, int pointid,Vec3D D );

    
    //  char str = new str[1000];
    char str1[256];
    char str2[256];
    
    FILE *fdtspoins;
    fdtspoins = fopen(file.c_str(), "r");
    char layer[256];
    float Lx,Ly,Lz;
    int NoPoints;
    
    if (fdtspoins == NULL){
        printf(" Error: Could not open file %s",file.c_str());
    }
    
    int readafile = fscanf(fdtspoins,"%s%s%s%d%s",str2,str2,str2,&NoPoints,str2);
    
    bool  check = fgets(str1, sizeof(str1), fdtspoins);
    check = fgets(str1, sizeof(str1), fdtspoins);

    
    std::vector<inclusion>  AllInc;
    
    
    float x,y,z;

    int id,tid,pid;
    for (int i=0;i<NoPoints;i++)
    {
        readafile = fscanf(fdtspoins,"%d%d%d%f%f%f",&id,&tid,&pid,&x,&y,&z);
        Vec3D D(x,y,z);
        inclusion p(id, tid,pid,D);
       // std::cout<<id<<"  "<<tid<<"  "<<pid<<"  "<<x<<"  "<<y<<"  "<<z<<"\n";
        AllInc.push_back(p);
    }
    
    
    return AllInc;
    
}
std::vector<exclusion> ReadDTSFolder::ReadExclusionObjects(std::string file)
{
    //    exclusion(int id, int pointid, double radius );
    
    
    //  char str = new str[1000];
    
    /// Read the header and find the number of the exclusion
    char str1[256];
    char str2[256];
    
    FILE *fdtspoins;
    fdtspoins = fopen(file.c_str(), "r");
    int NoPoints;
    
    if (fdtspoins == NULL){
        printf(" Error: Could not open file %s",file.c_str());
    }

    int readafile = fscanf(fdtspoins,"%s%s%s%d%s",str2,str2,str2,&NoPoints,str2);
    bool  check = fgets(str1, sizeof(str1), fdtspoins);
    check = fgets(str1, sizeof(str1), fdtspoins);
    ///
    
    std::vector<exclusion>  AllExc;
    
    
    float r;
    
    int id,pid;
    for (int i=0;i<NoPoints;i++)
    {
        readafile = fscanf(fdtspoins,"%d%d%f",&id,&pid,&r);
        exclusion p(id,pid,r);
        AllExc.push_back(p);
    }
    
    
    return AllExc;
    
}

bool ReadDTSFolder::FileExist (const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}
