

#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdlib>
#include "Curvature.h"
#include "In_OR_Out.h"
#include "WriteFiles.h"
#include "help.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "VMDOutput.h"
#include "Trajectory.h"
#include "Topology.h"
#include "Surface_Mosaicing.h"
#include "VertexMove.h"
#include "MakePBCTS.h"
#include "Traj_XXX.h"


/*
 this class is takes a TS file, ask for a point and checks if the point is inside the TS file.
 */
In_OR_Out::In_OR_Out(std::string filename)
{
    
    std::cout<<" Some fun, is a point inside or outside the TS \n";


    if(FileExist(filename) == false)
    {
        std::cout<<"Error: The TS file does not exist \n";
    }
    else
    {
        m_Box(0) = 3;
        m_Box(1) = 3;
        m_Box(2) = 3;
        m_pBox=&m_Box;
        Initialize(filename);   // read this TS file

        
        
        while(true)
        {
            std::cout<<" Enter 3 numbers for X, Y and Z coordinate of your point \n";
            double x,y,z;
            std::cin>>x>>y>>z;
            int inside = 0;
            for (std::vector<triangle *>::iterator it = m_pAllT.begin() ; it != m_pAllT.end(); ++it)
            {
                (*it)->UpdateNormal_Area(m_pBox);
                Vec3D N = (*it)->GetAreaVector();
                vertex *v = (*it)->GetV1();
                Vec3D X0(v->GetVXPos(),v->GetVYPos(),v->GetVZPos());
                Vec3D X(x,y,z);
                double value = N.dot((X-X0),N);
                if(value<0)
                {
                    inside--;
                }
                else if(value>0)
                {
                    inside++;
                }
                else
                {
                    std::cout<<" Point is on the surface \n";
                    std::cout<<0<<"\n";
                    break;
                }



            }
            if(inside>0)
            {
                std::cout<<" Point is outside \n";
                std::cout<<inside<<"\n";
            }
            else if(inside<0)
            {
                std::cout<<" Point is inside \n";
                std::cout<<inside<<"\n";
            }
            else
            {
                std::cout<<" unexpected \n";
            }
            
            
            std::cout<<" ****************** New try ************** \n";

        }


        

    }
}

In_OR_Out::~In_OR_Out()
{
    
}
void  In_OR_Out::UpdateGeometry()
{
    
    
    for (std::vector<triangle *>::iterator it = m_pAllT.begin() ; it != m_pAllT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    for (std::vector<links *>::iterator it = m_pHalfLinks1.begin() ; it != m_pHalfLinks1.end(); ++it)
    {
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    for (std::vector<vertex *>::iterator it = m_pAllV.begin() ; it != m_pAllV.end(); ++it)
        Curvature P(*it);
   
    
}
bool In_OR_Out::FileExist (const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}
void In_OR_Out::Initialize(std::string file)
{
    double lm = 0.1;
    double lmax =1000;
    double mina = 0;
    Topology S(m_pBox, &mina, &lm, &lmax);
    Traj_XXX TSI(m_pBox);


if(file.at(file.size()-1)=='q' && file.at(file.size()-2)=='.')
{

    S.FastReadQFile(file);
    bool topohealth= S.GetTopologyHealth();

    if(topohealth==false)
    {
        std::cout<<" error: Provided TS file is bad \n";
    }
    
    m_pAllV.clear();
    m_pInc.clear();
    m_pAllT.clear();
    m_pAllLinks.clear();
    m_pHalfLinks1.clear();
    m_pHalfLinks2.clear();
    m_pAllV=S.GetVertex();
    m_pAllT=S.GetTriangle();
    m_pAllLinks=S.GetLinks();
    m_pHalfLinks1=S.GetHalfLinks();
    m_pHalfLinks2=S.GetMHalfLinks();

    
}
else if(file.at(file.size()-1)=='i' && file.at(file.size()-2)=='s' && file.at(file.size()-3)=='t')
{
        m_pAllV.clear();
        m_pInc.clear();
        m_pAllT.clear();
        m_pAllLinks.clear();
        m_pHalfLinks1.clear();
        m_pHalfLinks2.clear();
        
        TSI.ReadTSI(file);
        m_pAllV=TSI.GetVertex();
        m_pAllT=TSI.GetTriangle();
        m_pAllLinks=TSI.GetLinks();
        m_pHalfLinks1=TSI.GetHalfLinks();
        m_pHalfLinks2=TSI.GetMHalfLinks();
        m_pInc=TSI.GetInclusion();
}
else
{
std::cout<<" error: Unknown TS File Format "<<file<<"\n";
}

    
    
    
}
