

#include <stdio.h>
#include "Wall.h"
#include "GenerateUnitCells.h"

Wall::Wall()
{

    m_State = false;
    m_Density = 1;
    m_H = 0;
    m_BeadName = "WL";
    m_CellSize = 3;
    m_Uniform = false;

}
Wall::~Wall()
{
    
}
void Wall::UpdateBox(Vec3D *x)
{
m_pBox=x;
}
void Wall::UpdateState(bool x)
{
    m_State = x;
}
void Wall::UpdateUniform(bool x)
{
    m_Uniform = x;
}
void Wall::UpdateCellSize(double x)
{
    m_CellSize = x;
}
void Wall::UpdateH(double x)
{
    m_H = x;
}
void Wall::UpdateDen(double x)
{
    m_Density = x;
}
void Wall::UpdateBeadName(std::string x)
{
    m_BeadName = x;
}
void Wall::PrintWallState()
{
    std::cout<<"---> Wall info \n";
    std::cout<<"---> state "<<m_State<<"   H   "<<m_H<<"   Density "<<m_Density<<"  Wall beads name "<<m_BeadName<<"\n";

 }
void Wall::CreateWall(std::vector<point*>  p1, std::vector<point*>  p2)
{
    if(m_State==true)
    {
        std::cout<<"We are creating a wall beads for this system with: \n";
        PrintWallState();
        Wall_itp();   //============ Make the wall itp file
        
        
        //===== We find largest possible wall density, which is equal to the 1/largest point area, we also create temprary beads for later use
        double maxarea = 0;
        int i=0;
        int resid = 0;
        for (std::vector<point*>::iterator it = p1.begin() ; it != p1.end(); ++it)
        {
                double area = (*it)->GetArea();
                if(maxarea<area)
                    maxarea = area;
        }
        i=0;
        resid = 0;
        for (std::vector<point*>::iterator it = p2.begin() ; it != p2.end(); ++it)
        {
            double area = (*it)->GetArea();
            if(maxarea<area)
                maxarea = area;
        }
        double maxden = 1.0/maxarea;
        m_Density = maxden*m_Density;
        
        std::cout<<"Note ----> Maximum wall density beads that can be created is "<<maxden<<" [particle]/nm^2, if you need more, use PLM and increase Mashno \n";
        //===== making wall beads of all availabe points

        std::vector<bead> b1 = MakeUniformBeads(p1);
        std::vector<bead> b2 = MakeUniformBeads(p2);
        

        
        
        for (std::vector<bead>::iterator it = b1.begin() ; it != b1.end(); ++it)
         {
              m_AllWallBeads.push_back(*it);
              m_AllWallPoints.push_back(p1.at((*it).GetID()));
         }
        
        for (std::vector<bead>::iterator it = b2.begin() ; it != b2.end(); ++it)
        {
            m_AllWallBeads.push_back(*it);
            m_AllWallPoints.push_back(p2.at((*it).GetID()));
        }

        std::cout<<"Note ----> number of the wall beads for the upper monolayer is: "<<b1.size()<<"\n";
        std::cout<<"Note ----> number of the wall beads for the inner monolayer is: "<<b2.size()<<"\n";

        

        
    }
    else
    {
        
    }
    
}
void Wall::Wall_itp()
{
    std::string file="Wall.itp";
    FILE *fitp;
    fitp = fopen(file.c_str(), "w");
    const char* atype=" [ atomtypes ] ";
    const char* AA=" A ";
    fprintf(fitp,  "%s\n",atype);
    fprintf(fitp, "%5s%8.3f%8.3f%5s%8.3f%8.3f\n",m_BeadName.c_str(),72.0,0.00,AA,0.0,0.0);
    
    {
        const char* tem=" [ nonbond_params ] ";
        fprintf(fitp,  "%s\n",tem);
        const char* tem2=" C1     WL     1    0.47  0.500 ";
        fprintf(fitp,  "%s\n",tem2);
    }
    {
        const char* tem=" [ moleculetype ] ";
        fprintf(fitp,  "%s\n",tem);
    }
    {
        const char* tem=" Wall              1 ";
        fprintf(fitp,  "%s\n",tem);
    }
    {
        const char* tem=" [ atoms ] ";
        fprintf(fitp,  "%s\n",tem);
    }
    {
        std::string str =  "1      "+m_BeadName+"      1        wall  "+m_BeadName +"     1       0";
        const char* tem=str.c_str();
        fprintf(fitp,  "%s\n",tem);
    }
    {
        const char* tem1="[ position_restraints ]";
        fprintf(fitp,  "%s\n",tem1);
        const char* tem2="1    1 1000 1000 1000";
        fprintf(fitp,  "%s\n",tem2);
    }
    fclose (fitp);
}
/// This function creats a uniform bead distribution based on the given density
std::vector<bead> Wall::MakeUniformBeads(std::vector<point*> &mypoints)
{


    std::vector<bead> ReturnB;
    int i=0;
    int resid = 0;
    std::vector<bead> TB;
    for (std::vector<point*>::iterator it = mypoints.begin() ; it != mypoints.end(); ++it)
    {
        
        double p_area = (*it)->GetArea();
      
        
        Vec3D X=(*it)->GetPos();
        /*if(m_Uniform==false)
        {
        Vec3D N=(*it)->GetNormal();
        X=X+N*(m_H);
        }*/
        bead  tb(i, m_BeadName, m_BeadName, "Wall", resid,X(0), X(1),X(2));
        TB.push_back(tb);
        i++;
        resid++;
        


    }

    std::vector<bead*> pTB;
    for (std::vector<bead>::iterator it = TB.begin() ; it != TB.end(); ++it)
        pTB.push_back(&(*it));
    
    
    
  if(m_Uniform==true)
  {
    GenerateUnitCells GN(pTB,m_pBox, 1, m_CellSize);
    GN.Generate();

    std::map <std::string, UnitCell > AUC = GN.GetAllCNTCells();
    

    for (std::map <std::string, UnitCell >::iterator it = AUC.begin(); it != AUC.end(); it++)
    {
        UnitCell tu= it->second;
        std::vector <bead *> B = tu.GetBeadList();
        double area = 0;
        for (std::vector<bead*>::iterator it = B.begin() ; it != B.end(); ++it)
        {
            area+=(mypoints.at((*it)->GetID()))->GetArea();
        }
        int newN =0;
        if(m_Density!=0 && B.size()!=0)
        newN=int(m_Density*area)+1;
      //  std::cout<<newN<<" area1 "<<m_Density<<"  "<<area<<" \n";

        if(newN>B.size())
        {
            newN =B.size();
            std::cout<<" warning, should not happen. Please report to developer with ID number wall3245643 \n";
        }
        double ratio = double(newN)/double(B.size());
        // === For developer, if the number of the beads are too little this will not give the exact number
        for (std::vector<bead*>::iterator it = B.begin() ; it != B.end(); ++it)
            (*it)->UpdateHasMol(false);
        
        

        int accepted = 0;
        while(newN>accepted)
        for (std::vector<bead*>::iterator it = B.begin() ; it != B.end(); ++it)
        {

            int trand = rand()%100000000;
            double prob=double(trand)/100000000.0;
            double p_area = (mypoints.at((*it)->GetID()))->GetArea();
          //  int id=i;
            Vec3D X=(mypoints.at((*it)->GetID()))->GetPos();
            Vec3D N=(mypoints.at((*it)->GetID()))->GetNormal();
            X=X+N*(m_H);
            (*it)->UpdatePos(X(0),X(1),X(2));
            if(prob<p_area/area && newN>accepted && (*it)->Hasmol()==false)
            {
                ReturnB.push_back(*(*it));
                accepted++;
                (*it)->UpdateHasMol(false);

            }
        }
        if (accepted<newN)
            std::cout<<newN-accepted<<"   "<<accepted<<" should have been "<<newN<<"\n";

        
     }
    }
    else
    {
        std::cout<<" here in the wall is ok \n";
        for (std::vector<bead>::iterator it = TB.begin() ; it != TB.end(); ++it)
        {

            ReturnB.push_back(*it);
            
        }
    }

    return ReturnB;
 
}











