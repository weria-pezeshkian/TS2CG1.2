  #if !defined(AFX_PointBasedBlueprint_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_PointBasedBlueprint_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "PointBasedBlueprint.h"
#include "Def.h"
#include "FlatPointMaker.h"
#include "Sphere.h"
#include "Cylinder.h"
#include "SHGeneric1DPBCPointMaker.h"

/*
 
 further extentions
 
 1). in random inc, we should prevent them to be close
 2.) angle and theta is not used yet
 3.) pattern based protein insertion

 */
PointBasedBlueprint::PointBasedBlueprint(Argument *pArgu)
{
    Nfunction f;      // In this class there are some useful function and we can use it.
    
    
    std::string function = pArgu->GetFunction();
    std::string dtsfoldername = pArgu->GetDTSFolder();              // get the name of the folder generated by the PLM
    std::string ifilename = pArgu->GetStructureFileName();          // get the input file name


    ReadDTSFolder DTSFolder;

    if(function=="analytical_shape")
    {
        std::string ftype= functiontype(ifilename);
    
        if(ftype == "Flat")
        {
            std::cout<<"---> note: Flat bilayer will be made \n";
            FlatPointMaker  Fu(pArgu);
            m_PointUp = Fu.GetUpPoint();
            m_PointDown = Fu.GetInPoint();
            m_WPointUp  = Fu.GetWallPoint1();
            m_WPointDown  = Fu.GetWallPoint2();
            m_Box=Fu.GetBox();
        }
         else if(ftype == "1D_PBC_Fourier")
        {
            std::cout<<"---> note: shape from 1D_PBC_Fourier will be made \n";
            SHGeneric1DPBCPointMaker  Fu(pArgu);
            m_PointUp = Fu.GetUpPoint();
            m_PointDown = Fu.GetInPoint();
            m_WPointUp  = Fu.GetWallPoint1();
            m_WPointDown  = Fu.GetWallPoint2();
            m_Box=Fu.GetBox();
        }
        else if(ftype == "Sphere")
        {
            std::cout<<"---> note: vesicle will be made \n";
            Sphere  Fu(pArgu);
            m_PointUp = Fu.GetUpPoint();
            m_PointDown = Fu.GetInPoint();
            m_WPointUp  = Fu.GetWallPoint1();
            m_WPointDown  = Fu.GetWallPoint2();
            m_Box=Fu.GetBox();
        }
        else if(ftype == "Cylinder")
        {
            std::cout<<"---> note: vesicle will be made \n";
            Cylinder  Fu(pArgu);
            m_PointUp = Fu.GetUpPoint();
            m_PointDown = Fu.GetInPoint();
            m_WPointUp  = Fu.GetWallPoint1();
            m_WPointDown  = Fu.GetWallPoint2();
            m_Box=Fu.GetBox();
        }
        else
        {
            std::cout<<"----> error: the shape defined in the str file is unknown :) \n";
            std::exit(0);
        }
    }
    else
    {
        
        DTSFolder.Read(dtsfoldername);
        m_PointUp = DTSFolder.GetUpperPoints();
        m_PointDown = DTSFolder.GetInnerPoints();
        m_Inc = DTSFolder.GetInclusion();
        m_Exc = DTSFolder.GetExclusion();
        m_Box= DTSFolder.GetBox();

    }
    
    if(m_PointDown.size()==0)
    m_monolayer = true;
    
    //==== puting the read data in a pointer container
        m_pBox =&m_Box;
     for (std::vector<point>::iterator it = m_PointUp.begin() ; it != m_PointUp.end(); ++it)
         m_pPointUp.push_back(&(*it));
    
    for (std::vector<point>::iterator it = m_PointDown.begin() ; it != m_PointDown.end(); ++it)
        m_pPointDown.push_back(&(*it));
    
    for (std::vector<inclusion>::iterator it = m_Inc.begin() ; it != m_Inc.end(); ++it)
        m_pInc.push_back(&(*it));
    
    for (std::vector<exclusion>::iterator it = m_Exc.begin() ; it != m_Exc.end(); ++it)
        m_pExc.push_back(&(*it));
    
    for (std::vector<inclusion*>::iterator it = m_pInc.begin() ; it != m_pInc.end(); ++it)
    {
        int pid=(*it)->GetPointID();
        
        if(pid<m_pPointUp.size() && pid>=0)
        {
            (m_pPointUp.at(pid))->UpdateInclusion(*it);
        }
        else
        {
            std::cout<<"--->error: the point id for the inclusion is out of range: \n ";
            std::cout<<" point id "<<pid<<" while # points is "<<m_pPointUp.size()<<"\n";
            std::exit(0);
        }
        
    }
    for (std::vector<point>::iterator it = m_WPointUp.begin() ; it != m_WPointUp.end(); ++it)
        m_pWPointUp.push_back(&(*it));
   for (std::vector<point>::iterator it = m_WPointDown.begin() ; it != m_WPointDown.end(); ++it)
       m_pWPointDown.push_back(&(*it));
    
//==== end puting the read data in a pointer container

    
    //=============== make wall; Wall info and data
        m_Wall = pArgu->GetWall();
        m_Wall.UpdateBox(m_pBox);
         if(function=="analytical_shape")
             m_Wall.CreateWall(m_pWPointUp,m_pWPointDown);
         else
             m_Wall.CreateWall(m_pPointUp,m_pPointDown);

        m_pWall = &m_Wall;

    //********* Lets exclude the points based on exclusion
    if(m_pExc.size()!=0)
    {
        std::cout<<" Note: we are excluding points based on exclusion, If it is slow, contact the developer \n";
        
        for ( std::vector<exclusion*>::iterator it = m_pExc.begin(); it != m_pExc.end(); it++ )
        {
            int pointid=(*it)->GetPointID();

            if(pointid<0 || pointid>m_pPointUp.size())
            {
                std::cout<<"---> error point id for the exclusion is wrong \n";
                std::exit(0);
            }

            double R = (*it)->GetRadius();
            point *Up_p1=m_pPointUp[pointid];
            Vec3D Pos = Up_p1->GetPos();
            Vec3D N = Up_p1->GetNormal();
            if(R!=0)
            Up_p1->UpdateArea(0);

            for ( std::vector<point*>::iterator it1 = m_pPointUp.begin(); it1 != m_pPointUp.end(); it1++ )
            {
                Vec3D Pos1 = (*it1)->GetPos();
                Vec3D DP = Pos1-Pos;
                double dist = DP.norm();
                Vec3D UnitDP =DP*(1/dist);
                double sinT = fabs((UnitDP*N).norm());
                double cosT = fabs(N.dot(UnitDP,N));

                if(dist*sinT<=R && dist*cosT<6)
                    (*it1)->UpdateArea(0);
  
            }
            for ( std::vector<point*>::iterator it1 = m_pPointDown.begin(); it1 != m_pPointDown.end(); it1++ )
            {
                Vec3D Pos1 = (*it1)->GetPos();
                Vec3D DP = Pos1-Pos;
                double dist = DP.norm();
                Vec3D UnitDP =DP*(1/dist);
                double sinT = fabs((UnitDP*N).norm());
                double cosT = fabs(N.dot(UnitDP,N));

                if(dist*sinT<=R && dist*cosT< 6)
                    (*it1)->UpdateArea(0);
                    
                
            }
        }
    }

}
PointBasedBlueprint::~PointBasedBlueprint()
{
    
}
std::string PointBasedBlueprint::functiontype(std::string filename)
{
    std::string ftype;
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
            else if(flag==true && Line.at(0)=="ShapeType")
            {

                    if(Line.size()<2)
                        std::cout<<" Error: ShapeType information in the str file is not correct \n";
                    else
                    ftype = Line.at(1);
                
                break;

            }
            
            
        }
        
        
    }
    
    
    file.close();
    

    return ftype;
}






#endif



