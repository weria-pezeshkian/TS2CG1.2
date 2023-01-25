  #if !defined(AFX_Membrane_Builder_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_Membrane_Builder_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "Membrane_Builder.h"
#include "GroFile.h"
#include "ReadDTSFolder.h"
#include "GenerateUnitCells.h"
#include "Def.h"
#include "PDBFile.h"
#include "GenDomains.h"
#include "FlatPointMaker.h"
#include "Sphere.h"
#include "Cylinder.h"



Membrane_Builder::Membrane_Builder(Argument *pArgu)
{
    std::cout<<"==============  PCG membrane builder: from function to membrane ======================"<<"\n";
    std::cout<<"****************** "<< SoftWareName <<" ******************   "<<"\n";
    std::cout<<"******************  Version:  "<<SoftWareVersion<<" ****************** \n";
    std::cout<<"=========================================================================================================="<<"\n";
    std::string ifilename = pArgu->GetStructureFileName();
    std::string ftype= functiontype(ifilename);
    Nfunction f;      // In this class there are some useful function and we can use it.
    if(ftype == "1D_PBC_Fourier")
    {
        std::cout<<"shape from 1D_PBC_Fourier will be made \n";
        SHGeneric1DPBCPointMaker  Fu(pArgu);
        m_Point1 = Fu.GetUpPoint();
        m_Point2 = Fu.GetInPoint();
        
        m_WallPoint1  = Fu.GetWallPoint1();
        m_WallPoint2  = Fu.GetWallPoint2();
        m_Box=Fu.GetBox();
        m_pBox =&m_Box;
    }
    else if(ftype == "Flat")
    {
        std::cout<<"Flat bilayer will be made \n";
        FlatPointMaker  Fu(pArgu);
        m_Point1 = Fu.GetUpPoint();
        m_Point2 = Fu.GetInPoint();
        m_WallPoint1  = Fu.GetWallPoint1();
        m_WallPoint2  = Fu.GetWallPoint2();
        m_Box=Fu.GetBox();
        m_pBox =&m_Box;
    }
    else if(ftype == "Sphere")
    {
        std::cout<<"vesicle will be made \n";
        Sphere  Fu(pArgu);
        m_Point1 = Fu.GetUpPoint();
        m_Point2 = Fu.GetInPoint();
        m_WallPoint1  = Fu.GetWallPoint1();
        m_WallPoint2  = Fu.GetWallPoint2();
        m_Box=Fu.GetBox();
        m_pBox =&m_Box;
    }
    else if(ftype == "Cylinder")
    {
        std::cout<<"vesicle will be made \n";
        Cylinder  Fu(pArgu);
        m_Point1 = Fu.GetUpPoint();
        m_Point2 = Fu.GetInPoint();
        m_WallPoint1  = Fu.GetWallPoint1();
        m_WallPoint2  = Fu.GetWallPoint2();
        m_Box=Fu.GetBox();
        m_pBox =&m_Box;
    }
    else
    {
        std::cout<<"-> aborted! ---> the shape defined in the str file is unknown :) \n";
        std::exit(0);
    }
    m_monolayer = false;
    srand (pArgu->GetSeed());
    //===================== OutPut file name declaration and finding input file names ========================================
    std::string gname = pArgu->GetGeneralOutputFilename();
    m_FinalOutputGroFileName =gname+".gro";
    std::string m_FinalTopologyFileName=gname+".top";
    std::string strfilename = pArgu->GetStructureFileName();
    m_ResID = 1;
    m_Renormalizedlipidratio = pArgu->GetRenorm();
    m_Iter = pArgu->GetIter();
    GenerateMolType  MOLTYPE(pArgu);
    m_MoleculesType = MOLTYPE.GetMolType();
    if(MOLTYPE.GetHealth()==false)
    {
        std::cout<<"-> aborted! You are allowed to try one more time. Kidding, please do not :) \n";
        std::exit(0);
        
    }
    
    
    //********************** Finding the total area of the layers
    for (std::vector<point >::iterator it = m_Point1.begin() ; it != m_Point1.end(); ++it)
    {
        m_pPoint1.push_back(&(*it));
        m_TotalAreaUp+=it->GetArea();
    }
    for (std::vector<point >::iterator it = m_Point2.begin() ; it != m_Point2.end(); ++it)
    {
        m_pPoint2.push_back(&(*it));
        m_TotalAreaDown+=it->GetArea();
    }
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
    CWall.UpdateBox(m_pBox);
    for (std::vector<point >::iterator it = m_WallPoint1.begin() ; it != m_WallPoint1.end(); ++it)
    {
        m_pWallPoint1.push_back(&(*it));
    }
    for (std::vector<point >::iterator it = m_WallPoint2.begin() ; it != m_WallPoint2.end(); ++it)
    {
        m_pWallPoint2.push_back(&(*it));
    }
    CWall.CreateWall(m_pWallPoint1,m_pWallPoint2);
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
Membrane_Builder::~Membrane_Builder()
{
    
}
void Membrane_Builder::GenLipid(MolType moltype, int listid, Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2)
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
double Membrane_Builder::dist2between2Points(Vec3D X1,Vec3D X2)
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
void Membrane_Builder::WriteFinalGroFile()
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

Tensor2 Membrane_Builder::Rz(double cos, double sin)
{
    Tensor2  R;
    
    R(0,0) = cos;
    R(0,1) = -sin;
    R(1,0) = sin;
    R(1,1) = cos;
    R(2,2) = 1;

    
    
    return R;
}
/*void Membrane_Builder::CreateWallBead(std::vector<point*>  p1, std::vector<point*>  p2)
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
Tensor2  Membrane_Builder::TransferMatLG(Vec3D Normal, Vec3D t1, Vec3D t2)
{
    
    Tensor2  GL(t1,t2,Normal);
    Tensor2 LG=GL.Transpose(GL);
    return LG;
    
}
std::string Membrane_Builder::functiontype(std::string filename)
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
        
        std::cout<<str<<"\n";
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



