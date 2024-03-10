  #if !defined(AFX_Shape_1DSinBuilder_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_Shape_1DSinBuilder_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "Shape_1DSinBuilder.h"
#include "GroFile.h"
#include "ReadDTSFolder.h"
#include "GenerateUnitCells.h"
#include "Def.h"
#include "PDBFile.h"
#include "GenDomains.h"
Shape_1DSinBuilder::Shape_1DSinBuilder(Argument *pArgu)
{
    
    m_State = pArgu->Get1DSinState();
    m_monolayer = false;
    
    srand (pArgu->GetSeed());
    std::cout<<"\n";
    std::cout<<"=========================================================================================================="<<"\n";
    std::cout<<"****************** "<< SoftWareName <<" ******************   "<<"\n";
    std::cout<<"******************  Version:  "<<SoftWareVersion<<" ****************** \n";
    std::cout<<"=========================================================================================================="<<"\n";

    Nfunction f;      // In this class there are some useful function and we can use it.


    
    //===================== OutPut file name declaration and finding input file names ========================================
    std::string gname = pArgu->GetGeneralOutputFilename();
    m_FinalOutputGroFileName =gname+".gro";
    std::string m_FinalTopologyFileName=gname+".top";
    std::string strfilename = pArgu->GetStructureFileName();
    m_ResID = 1;
    m_Renormalizedlipidratio = pArgu->GetRenorm();
    
    m_APLLipids = m_State.APL ;
    m_APLWall  = m_State.APW;
    
    m_Iter = pArgu->GetIter();
    //==========================================================================================================
    GenerateMolType  MOLTYPE(pArgu);
    m_MoleculesType = MOLTYPE.GetMolType();
    
    if(MOLTYPE.GetHealth()==false)
    {
        std::cout<<"-> aborted! You are allowed to try one more time. Kidding, please do not :) \n";
        std::exit(0);
        
    }
    m_Box(0)=m_State.Lx;     m_Box(1)=m_State.Ly;     m_Box(2)=m_State.Lz;
    m_pBox = (&m_Box);
    //********************** Finding the total area of the layers

    m_TotalAreaUp = CalculateArea_MakePoints(1,m_APLLipids);
    m_TotalAreaDown = CalculateArea_MakePoints(-1,m_APLLipids);

    for (int i=0;i<m_Point1.size();i++)
        m_pPoint1.push_back(&(m_Point1.at(i)));
    
    for (int i=0;i<m_Point2.size();i++)
        m_pPoint2.push_back(&(m_Point2.at(i)));
    

    //TEST();
    if(m_pPoint2.size()==0)
        m_monolayer = true;

    if(m_monolayer == false)
    std::cout<<"--> Note: the total upper monolayer area is "<<m_TotalAreaUp<<" and the total lower monolayer area is "<<m_TotalAreaDown<<"\n";
    else
    std::cout<<"--> Note: the total monolayer area is "<<m_TotalAreaUp<<" \n";

    
    
    std::cout<<" now is time to add the lipids \n";
    // Make all the domain containing different lipids
    GenDomains GENDOMAIN(strfilename,m_pPoint1,m_pPoint2,m_Renormalizedlipidratio);
    std::vector<Domain*> pAllDomain = GENDOMAIN.GetDomains();
    std::cout<<" Number of the domains defined in the input file  "<< pAllDomain.size()/2 <<"\n";
    int layer = 0;
        std::cout<<"***************************  we aim to generate  ********************** \n";
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {

        layer++;
        std::vector<DomainLipid> DL = (*it)->GetDomainLipids();
        
        if(layer%2!=0 || m_monolayer == false)
        {
        std::cout <<"*     For domain with ID "<<(*it)->GetDomainID() <<" \n";
        for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
        {
            std::cout <<"*     "<<(*it2).MaxNo<<"  "<<(*it2).Name<<"     "<<std::endl ;

        }
        }
        if(layer%2!=0 && m_monolayer == false)
        std::cout <<"   in the upper monolayer \n";
        else if(layer%2==0 && m_monolayer == false)
        std::cout <<"   in the lower monolayer \n";
        if(layer%2!=0 && m_monolayer == true)
        std::cout <<"   in the  monolayer \n";
    }
    

    
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {
            std::vector<DomainLipid*> DL = (*it)->GetpDomainLipids();
       /*     std::vector<DomainLipid> DLL = (*it)->GetDomainLipids();

            int s=0;
            for ( std::vector<DomainLipid*>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
            {
                std::cout <<s<<" *   begining   "<<(*it2)->MaxNo<<"  "<<(*it2)->Name<<"     "<<std::endl ;
                s++;
            }
         s=0;
        for ( std::vector<DomainLipid>::iterator it2 = DLL.begin(); it2 != DLL.end(); it2++ )
        {
            std::cout <<s<<" *   beginingQ   "<<(*it2).MaxNo<<"  "<<(*it2).Name<<"     "<<std::endl ;
            s++;
        }*/
        int iteration = 0;
        int madetotallipids = 0;
        int lipidlistID = 0;
        int NoMadeLipid = 0 ;  // temperory because our strcture cannot increase it
        int CreateTotalLipids = (*it)->GetDomainTotalLipid();

        std::vector<point*>  dpoint = (*it)->GetDomainPoint();
        while(true)
        {
            iteration++;

                if( lipidlistID>DL.size()-1)
                    break;
                if(madetotallipids==CreateTotalLipids)
                    break;
                if(iteration>m_Iter*(dpoint.size()))
                {
                std::cout<<" Warning: With "<< m_Iter <<" iterations, we could not place the expected number of the lipids \n";
                std::cout<<" if you are unhappy, increase the number of the iteration with option -iter, or regenerate the points \n";
                break;
                }



            
            int pointid = rand()%(dpoint.size());
            point* tempoint = dpoint.at(pointid);

            DomainLipid *LL=DL.at(lipidlistID);
            int RNG=(rand()%CreateTotalLipids)+1;
            Vec3D  Dir(0,0,0);
            Vec3D N = tempoint->GetNormal();
            Vec3D T1 =   tempoint->GetP1();
            Vec3D T2 =   tempoint->GetP2();

            std::string ltype = LL->Name;
            Vec3D Pos = tempoint->GetPos();
            double area = tempoint->GetArea();
            double rn = double(rand()%(1000000))/1000000.0;
            double prob=area/(LL->Ap);
           // std::cout<<DL.size()<<"  "<<LL->Ap<<"We 2222get here \n";

            if(prob>rn && LL->MaxNo>NoMadeLipid)
            {
                madetotallipids++;
                NoMadeLipid++;
                int t = LL->no_created;
                ((*it)->GetpDomainLipids()).at(lipidlistID)->no_created=t+1;
                //std::cout<< ((*it)->GetpDomainLipids()).at(lipidlistID)->no_created<<" "<<t+1<<"no beads \n";
                if (m_MoleculesType.count(ltype) == 0)
                    std::cout << "Error:-----> molecule name " <<ltype<<" does not exist in the lib files \n";
                
                GenLipid(m_MoleculesType.at(ltype), 0, Pos, N, Dir, T1, T2);
                (dpoint.at(pointid))->UpdateArea(0);
            }
            if(LL->MaxNo==NoMadeLipid )
            {
           //  std::cout<<"domain id "<<(*it)->GetDomainID()<<"  name "<<LL->Name<<" created for the domain "<<madetotallipids<<"  domain "<<CreateTotalLipids<<"  "<<lipidlistID<<"  "<<NoMadeLipid<<"  "<<LL->MaxNo<<"\n";
                lipidlistID++;
                NoMadeLipid = 0;


            }
            
        }

        
    }
    
    std::cout<<"*************************** Number of created Lipids,   ********************** \n";
    layer = 0;
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {
        layer++;
        std::vector<DomainLipid> DL = (*it)->GetDomainLipids();
        
        if(layer%2!=0 || m_monolayer == false)
        {
            std::cout <<"*     For domain with ID "<<(*it)->GetDomainID() <<" \n";
            for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
            {
                std::cout <<"*     "<<(*it2).no_created<<"  "<<(*it2).Name<<"     "<<std::endl ;
                
            }
        }
        if(layer%2!=0 && m_monolayer == false)
        std::cout <<"   in the upper monolayer \n";
        else if(layer%2==0 && m_monolayer == false)
        std::cout <<"   in the lower monolayer \n";
        if(layer%2!=0 && m_monolayer == true)
        std::cout <<"   in the  monolayer \n";
    }

    //=============== make wall; Wall info and data
    Wall CWall = pArgu->GetWall();
    m_Point1.clear();
    m_pPoint1.clear();
    m_pPoint2.clear();
    m_Point2.clear();
    m_State.H = m_State.H+CWall.GetH();

    
    CalculateArea_MakePoints(1,m_APLWall);
    CalculateArea_MakePoints(-1,m_APLWall);
    
    for (int i=0;i<m_Point1.size();i++)
        m_pPoint1.push_back(&(m_Point1.at(i)));
    
    for (int i=0;i<m_Point2.size();i++)
        m_pPoint2.push_back(&(m_Point2.at(i)));
    
    CWall.UpdateBox(m_pBox);
    CWall.CreateWall(m_pPoint1,m_pPoint2);
    std::vector<bead> WB = CWall.GetWallBead();
    std::vector<point*> WPoint = CWall.GetWallPoint();
    
    
    
    //==========  We write the wall at the end
    //=============== write the wall info
    if(WPoint.size()>0 && CWall.GetState()==true)
    {
        PDBFile pdb;
        std::string pdbfile = "wall.pdb";
        pdb.WritePDBFile(pdbfile, WPoint);
    }
    else if(CWall.GetState()==true)
    {
        std::cout<<"Note ----> No wall.pdb file will be generated since the total created wall beads are zero \n";
    }
    for (std::vector<bead>::iterator it = WB.begin() ; it != WB.end(); ++it)
    {
        m_FinalBeads.push_back((*it));
    }
    //============== End Wall info and data

    
    WriteFinalGroFile();
    
    
    //==========================================================================================================
    //=============== Open Topology files and make mols
    //==========================================================================================================
    std::ofstream Topgro;
    Topgro.open(m_FinalTopologyFileName.c_str());
    Topgro<<" ;This file was generated by TS Back Mapping \n";
    Topgro<<" [ system ] \n";
    Topgro<<" Expect a large membrane \n";
    Topgro<<" [ molecules ] \n";
    


    layer = 0;
    
    
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {
        layer++;
        std::vector<DomainLipid> DL = (*it)->GetDomainLipids();
        
        if(layer%2!=0 || m_monolayer == false)
        {
            Topgro <<"; domain "<<(*it)->GetDomainID() <<" \n";
            if(layer%2!=0 && m_monolayer == false)
            Topgro  <<" ;  in the upper monolayer \n";
            else if(layer%2==0 && m_monolayer == false)
            Topgro <<" ;  in the lower monolayer \n";
            for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
            {
                Topgro <<"     "<<(*it2).Name<<"  "<<(*it2).no_created<<"     "<<std::endl ;
                
            }
        }

    }
    if(WB.size()!=0)
    Topgro<<"Wall    "<<WB.size()<<"\n";



}
Shape_1DSinBuilder::~Shape_1DSinBuilder()
{
    
}
void Shape_1DSinBuilder::GenLipid(MolType moltype, int listid, Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2)
{
    //

    
    
     Tensor2 LG = TransferMatLG(Normal, t1, t2);
    
    
    
      std::vector<bead> vbeads = moltype.Beads;
            for ( std::vector<bead>::iterator it = vbeads.begin(); it != vbeads.end(); it++ )
            {
                Vec3D BPos((*it).GetXPos(),(*it).GetYPos(),(*it).GetZPos());
                Vec3D vX = LG*BPos+ Pos;
                int beadid = (m_FinalBeads.size()+1);
                bead TemB(beadid, (*it).GetBeadName(), (*it).GetBeadType(), (*it).GetResName(), m_ResID, vX(0), vX(1),vX(2));
                m_FinalBeads.push_back(TemB);
            }
    
    m_ResID++;
    
    



    
}
double Shape_1DSinBuilder::dist2between2Points(Vec3D X1,Vec3D X2)
{
    
    double dist2=0;
    
    double x1=X1(0);
    double y1=X1(1);
    double z1=X1(2);
    
    double x2=X2(0);
    double y2=X2(1);
    double z2=X2(2);
    
    
    double dx=x2-x1;
    double dy=y2-y1;
    double dz=z2-z1;
    
    if(fabs(dx)>(*m_pBox)(0)/2.0)
    {
        if(dx<0)
            dx=(*m_pBox)(0)+dx;
        else if(dx>0)
            dx=dx-(*m_pBox)(0);
    }
    if(fabs(dy)>(*m_pBox)(1)/2.0)
    {
        if(dy<0)
            dy=(*m_pBox)(1)+dy;
        else if(dy>0)
            dy=dy-(*m_pBox)(1);
    }
    if(fabs(dz)>(*m_pBox)(2)/2.0)
    {
        if(dz<0)
            dz=(*m_pBox)(2)+dz;
        else if(dz>0)
            dz=dz-(*m_pBox)(2);
    }

    dist2=dx*dx+dy*dy+dz*dz;
    return dist2;
}
void Shape_1DSinBuilder::WriteFinalGroFile()
{
    
   

    FILE *fgro;
    fgro = fopen(m_FinalOutputGroFileName.c_str(), "w");
    
    
    /// resid  res name   noatom   x   y   z
    const char* Title=" System ";
    int Size=m_FinalBeads.size();
    
    fprintf(fgro,  "%s\n",Title);
    fprintf(fgro, "%5d\n",Size);
    int i=0;
    for (std::vector<bead>::iterator it = m_FinalBeads.begin() ; it != m_FinalBeads.end(); ++it)
    {
        
        i++;
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        
        
        const char* A1=((*it).GetResName()).c_str();
        const char* A2=((*it).GetBeadName()).c_str();
        int resid=(*it).GetResid();
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid%100000,A1,A2,i%100000,x,y,z );
        
    }
    
    
    fprintf(fgro,  "%10.5f%10.5f%10.5f\n",m_Box(0),m_Box(1),m_Box(2) );
    fclose(fgro);
    
    
}

Tensor2 Shape_1DSinBuilder::Rz(double cos, double sin)
{
    Tensor2  R;
    
    R(0,0) = cos;
    R(0,1) = -sin;
    R(1,0) = sin;
    R(1,1) = cos;
    R(2,2) = 1;

    
    
    return R;
}
/*void Shape_1DSinBuilder::CreateWallBead(std::vector<point*>  p1, std::vector<point*>  p2)
{
    std::string file="Wall.gro";
    FILE *fgro;
    fgro = fopen(file.c_str(), "w");
    
    
    /// resid  res name   noatom   x   y   z
    const char* Title=" Wall of Water Beads ";
    int Size=p1.size()+p2.size();
    
    fprintf(fgro,  "%s\n",Title);
    fprintf(fgro, "%5d\n",Size);
    int i=0;
    for (std::vector<point*>::iterator it = p1.begin() ; it != p1.end(); ++it)
    {
        
        i++;
        Vec3D X=(*it)->GetPos();
        Vec3D N=(*it)->GetNormal();
        
        X=X+N*(0.7);
        const char* A1="W2";
        const char* A2="PO4";
        int resid=1;
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid%10000,A1,A2,i%10000,X(0),X(1),X(2) );
        
    }
    for (std::vector<point*>::iterator it = p2.begin() ; it != p2.end(); ++it)
    {
        
        i++;
        Vec3D X=(*it)->GetPos();
        Vec3D N=(*it)->GetNormal();
        
        X=X+N*(0.7);
        const char* A1="W2";
        const char* A2="PO4";
        int resid=1;
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid%10000,A1,A2,i%10000,X(0),X(1),X(2) );
        
    }
    
    fprintf(fgro,  "%10.5f%10.5f%10.5f\n",m_Box(0),m_Box(1),m_Box(2) );
    fclose(fgro);
    
}
 */
Tensor2  Shape_1DSinBuilder::TransferMatLG(Vec3D Normal, Vec3D t1, Vec3D t2)
{
    
    Tensor2  GL(t1,t2,Normal);
    Tensor2 LG=GL.Transpose(GL);
    return LG;
    
}
double Shape_1DSinBuilder::CalculateArea_MakePoints(int layer, double APL)
{
    int N = m_State.Omega;
    double A = m_State.A;
    double H = m_State.H;
    double area = 0;
    double pi = acos(-1);
    double omega = double(N)*2*pi/((*m_pBox)(0));
    double dt = 0.0001/(*m_pBox)(1);
    int Nup = int((*m_pBox)(0)/dt+H/(2*dt))+1;
    int Ndown = int(H/(2*dt))+1;
    std::vector<Vec3D> Pos,NormalV,T1vec,T2vec;
    std::vector<Vec3D> V1,V2,V3,V4;

    Vec3D t2(0,1,0);
    std::vector <double> c;
    c.push_back(0);
    c.push_back(0);
    for (int i=-Ndown;i<Nup;i++)
    {
        double t = i*dt;
        Vec3D X = F(t,layer);

        
        if(X(0)>=0 && X(0)<(*m_pBox)(0))
        {
            V1.push_back(X);  //             Pos.push_back(X);
            V2.push_back(Normal(t,layer));  // NormalV.push_back(Normal(t,layer));
            V3.push_back(T1(t,layer));
            V4.push_back(t2);

        }

    }
    for (int i=1;i<V1.size();i++)
    {
        double DX = (V1.at(i-1))(0)-(V1.at(i))(0);
        double DZ = (V1.at(i-1))(2)-(V1.at(i))(2);
        area+=sqrt(DX*DX+DZ*DZ);
    }
    double Lenght = area;
    area = area*(*m_pBox)(1);
////////////////////////
////////
///////////////////////
    double DL = sqrt(APL);


    DL = Lenght/double(int(Lenght/DL));
    Vec3D DXOLD = V1.at(0);
    Pos.push_back(V1.at(0));
    NormalV.push_back(V2.at(0));
    T1vec.push_back(V3.at(0));
    T2vec.push_back(V4.at(0));
    for(int i=1;i<V1.size();i++)
   {
       Vec3D DX = DXOLD-V1.at(i);
       if(DX.norm()>=DL)
       {
           DXOLD = V1.at(i);
           
           Pos.push_back(V1.at(i));
           NormalV.push_back(V2.at(i));
           T1vec.push_back(V3.at(i));
           T2vec.push_back(V4.at(i));
       }

        
    }
    
    
    
    
    
/////////////////////////
    int beadid = 0;
    double dy = sqrt(APL);
    int NY = (*m_pBox)(1)/dy;
  //  dy = (*m_pBox)(1)/double(NY); this is good for even 

    APL = area/double(Pos.size()*NY);
    for (int i=0;i<Pos.size();i++)
    {
        for (int j=0;j<NY;j++)
        {
            double y=(double(j))*dy+0.5*double(i%2)*dy;
            (Pos.at(i))(1) = y;
            point p(beadid, APL, Pos.at(i), NormalV.at(i), T1vec.at(i), T2vec.at(i) , c );
            if(layer==1)
                m_Point1.push_back(p);
            else if(layer==-1)
                m_Point2.push_back(p);
            beadid++;
        }
    }
    
  //  std::cout<<Pos.size()<<" APL  "<<area/double(Pos.size()*NY)<<"  "<<APL<<"\n";

    
    return area;
}
Vec3D Shape_1DSinBuilder::F(double t, int layer)
{
    int N = m_State.Omega;
    double A = m_State.A;
    double H = m_State.H;
    double area = 0;
    double pi = acos(-1);
    double omega = double(N)*2*pi/((*m_pBox)(0));
    
    Vec3D vec;
    
    vec(0) = t-layer*0.5*H*A*omega*cos(omega*t)/sqrt(1+A*omega*A*omega*cos(omega*t)*cos(omega*t));
    vec(2) = (*m_pBox)(2)/2+A*sin(omega*t)+layer*0.5*H/sqrt(1+A*omega*A*omega*cos(omega*t)*cos(omega*t));

    
    
    return vec;
}
Vec3D Shape_1DSinBuilder::Normal(double t, int layer)
{
    int N = m_State.Omega;
    double A = m_State.A;
    double H = m_State.H;
    double area = 0;
    double pi = acos(-1);
    double omega = double(N)*2*pi/((*m_pBox)(0));
    
    Vec3D vec;
    
    vec(0) = -layer*A*omega*cos(omega*t)/sqrt(1+A*omega*A*omega*cos(omega*t)*cos(omega*t));
    vec(2) = layer/sqrt(1+A*omega*A*omega*cos(omega*t)*cos(omega*t));
    
    
    
    return vec;
}
Vec3D Shape_1DSinBuilder::T1(double t, int layer)
{
    int N = m_State.Omega;
    double A = m_State.A;
    double H = m_State.H;
    double area = 0;
    double pi = acos(-1);
    double omega = double(N)*2*pi/((*m_pBox)(0));
    
    Vec3D vec;
    
    vec(0) = -A*omega*cos(omega*t)/sqrt(1+A*omega*A*omega*cos(omega*t)*cos(omega*t));
    vec(1) = -1/sqrt(1+A*omega*A*omega*cos(omega*t)*cos(omega*t));
    
    
    
    return vec;
}






#endif



