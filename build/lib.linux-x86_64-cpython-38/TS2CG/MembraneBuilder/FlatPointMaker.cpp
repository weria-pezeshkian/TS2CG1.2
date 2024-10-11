  #if !defined(AFX_FlatPointMaker_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_FlatPointMaker_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "FlatPointMaker.h"
#include "GroFile.h"
#include "ReadDTSFolder.h"
#include "GenerateUnitCells.h"
#include "Def.h"
#include "PDBFile.h"
#include "GenDomains.h"
FlatPointMaker::FlatPointMaker(Argument *pArgu)
{
    m_WallBox.push_back(0);
    m_WallBox.push_back(1);
    m_WallBox.push_back(0);
    m_WallBox.push_back(1);

    m_Box(0) = 10;m_Box(1) = 10;m_Box(2) = 10;
    m_pBox = &m_Box;
    m_Thickness = 4;
    m_Density = 2;
    m_WallDensity = 1;
    std::string ifilename = pArgu->GetStructureFileName();
    Initialize(ifilename);


    //=============== We need the H wall distance, and the only info needed about the wall here
    Wall CWall = pArgu->GetWall();
    double Hwall = CWall.GetH();

    //*********

    m_Point1 = CalculateArea_MakePoints(1, 1/m_Density,m_Thickness/2,false);
    m_Point2 = CalculateArea_MakePoints(-1, 1/m_Density,m_Thickness/2,false);

    
    m_WallPoint1 = CalculateArea_MakePoints(1, 1/m_WallDensity,m_Thickness/2+Hwall,true);
    m_WallPoint2 = CalculateArea_MakePoints(-1, 1/m_WallDensity,m_Thickness/2+Hwall,true);

}
FlatPointMaker::~FlatPointMaker()
{
    
}
std::vector<point> FlatPointMaker::CalculateArea_MakePoints(int layer, double APL,double H, bool wall)
{
    std::vector<point> Cpoints;
    std::vector <double> c;  // we set curvature to zero as it is not important for this function
    c.push_back(0);
    c.push_back(0);
    std::vector<Vec3D> Pos,NormalV,T1vec,T2vec;
    std::vector<Vec3D> V1, V2,V3,V4;
    double xmin = 0;
    double ymin = 0;
    double xmax = (*m_pBox)(0);
    double ymax = (*m_pBox)(1);
    
    if(wall==true)
    {
        xmin = (m_WallBox.at(0))*((*m_pBox)(0));
        xmax = (m_WallBox.at(1))*((*m_pBox)(0));
        ymin = (m_WallBox.at(2))*((*m_pBox)(1));
        ymax = (m_WallBox.at(3))*((*m_pBox)(1));
        
        
    }
    double area = (*m_pBox)(0)*(*m_pBox)(1);
    double Lenght = (*m_pBox)(0);
    double DL = sqrt(APL);
    double Nx = int((*m_pBox)(0)/DL)+1;
    double Ny = int((*m_pBox)(1)/DL)+1;
    double dx = (*m_pBox)(0)/double(Nx);
    double dy = (*m_pBox)(1)/double(Ny);

    for (int i=0;i<Nx;i++)
    {
        double x = dx*double(i);
        for (int j=0;j<Ny;j++)
        {
            double y = dy*double(j);

            Vec3D X (x,y,((*m_pBox)(2)/2)+layer*H);
            V1.push_back(X);
            Vec3D T1(double(1+layer)/2,double(1-layer)/2,0);
            Vec3D T2(double(1-layer)/2,double(1+layer)/2,0);
            Vec3D N(0,0,layer);
            V3.push_back(T1);
            V4.push_back(T2);
            V2.push_back(N);

        }
        
        
    }
    for(int i=0;i<V1.size();i++)
    {
            Pos.push_back(V1.at(i));
            NormalV.push_back(V2.at(i));
            T1vec.push_back(V3.at(i));
            T2vec.push_back(V4.at(i));
    }
    int beadid = 0;
    for (int i=0;i<Pos.size();i++)
    {
            if((Pos.at(i))(0)>=xmin && (Pos.at(i))(0)<=xmax && (Pos.at(i))(1)>=ymin && (Pos.at(i))(1)<=ymax)
            {
                point p(beadid, APL, Pos.at(i), NormalV.at(i), T1vec.at(i), T2vec.at(i) , c );
                Cpoints.push_back(p);
                beadid++;
            }
    }
    for (std::vector<point >::iterator it = Cpoints.begin() ; it != Cpoints.end(); ++it)
    {
        it->UpdateArea(area/double(Cpoints.size()));
    }
    
    
    return Cpoints;
    
}
void FlatPointMaker::Initialize(std::string filename)
{
    
    bool OK=true;
    Nfunction f;
    std::ifstream file;
    file.open(filename.c_str());
    bool flag = false;
    std::string str;

    while (true)
    {
        std::getline (file,str);
        if(file.eof())
            break;
        
        std::vector<std::string> Line = f.split(str);
        if(Line.size()!=0 && (Line.at(0)).at(0)!=';')
        {
            if((Line.at(0)).at(0)=='[' && flag==false)
            {
                str = f.trim(str);
                str.erase(std::remove(str.begin(), str.end(), '['), str.end());
                str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
                str.erase(std::remove(str.begin(), str.end(), ' '), str.end());

                if(str=="ShapeData")
                    flag = true;
            }
            else if((Line.at(0))=="End" && flag==true)
            {
                flag=false;
            }
            else if(flag==true)
            {
                if(Line.at(0)=="Box")
                {
                    if(Line.size()<4)
                        std::cout<<" Error: Box information in the str file is not correct \n";
                    else
                    {
                        m_Box(0) = f.String_to_Double(Line.at(1));
                        m_Box(1) = f.String_to_Double(Line.at(2));
                        m_Box(2) = f.String_to_Double(Line.at(3));

                    }
                }
                else if(Line.at(0)=="ShapeType")
                {
                }
                else if(Line.at(0)=="WallRange")
                {
                    if(Line.size()<5)
                        std::cout<<" Error: Wall range information in the str file is not correct \n";
                    else
                    {
                        m_WallBox.at(0) = f.String_to_Double(Line.at(1));
                        m_WallBox.at(1) = f.String_to_Double(Line.at(2));
                        m_WallBox.at(2) = f.String_to_Double(Line.at(3));
                        m_WallBox.at(3) = f.String_to_Double(Line.at(4));

                    }
                }
                else if(Line.at(0)=="Density")
                {
                    if(Line.size()<3)
                        std::cout<<" Error: Density information in the str file is not correct \n";
                    else
                    {
                        m_Density = f.String_to_Double(Line.at(1));
                        m_WallDensity = f.String_to_Double(Line.at(2));
                        
                    }
                }
                else if(Line.at(0)=="Thickness")
                {
                    if(Line.size()<2)
                        std::cout<<" Error: Thickness information in the str file is not correct \n";
                    else
                    {
                        m_Thickness = f.String_to_Double(Line.at(1));
                        
                    }
                }
            }


        }
            
            
        }
        
    
    file.close();

    
}






#endif



