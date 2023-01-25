 #if !defined(AFX_Sin_1D_CPP_7F4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_Sin_1D_CPP_7F4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "Sin_1D.h"
#include "GroFile.h"
#include "IndexFileReader.h"
#include  "IndexGroup.h"
#include  "GenerateUnitCells.h"
#include "PDBFile.h"


Sin_1D::Sin_1D(Argument *pArgu)
{


//////////  Atom number in the input file is the number of the atom in the coordinate file and it start from 1.... but index which is used in reading xtc started from 0
///////// because that we use Atomno-1



std::cout<<"\n";
std::cout<<"------------- Sin_1D is running: Which mean your inputs were ok, but let see if its make sense too  -------------"<<"\n";
    std::cout<<"------------- Please note: time unit is ns in this software  -------------"<<"\n";

std::cout<<"\n";


Nfunction f;      // In this class there are some useful function and we can use it.

    /// read data from argument class
    
    
    std::string pdbfilename = pArgu->GetpdbFileName();
    PDBFile pdb;
    m_WallPoints = pdb.ReadPDBFile(pdbfilename);

    std::string trjfilename = pArgu->GetTrjFile();
    int inialtime = pArgu->GetBStep();
    int finaltime = pArgu->GetEStep();
    std::string gfilename = pArgu->GetGeneralOutputFilename();
    int framejump = pArgu->GetJump();
    std::string indexfilename = pArgu->GetIndexFileName();
    double cuttoff=pArgu->GetCutoff();
    int NoFrame = 0;
    int acceptedframe = 0;
    bool isfirstframe = true;
    int frame0time = 0;
    std::string grofilename = pArgu->GetGroFileName();
    GroFile Gro(grofilename);
    m_pBox = Gro.GetBox();
    m_pAllBeads = Gro.GetpAllBeads();
    std::cout<<"  System contains "<<m_pAllBeads.size()<<"   beads \n";
    GenerateUnitCells   GUN(m_pAllBeads,pArgu,m_pBox);
    
    std::ofstream timeseri;
    timeseri.open("meanC_Corr.xvg");
    
    timeseri<<"@@     time (ns)    average_Mean_Curvature  variance "<<"\n";
    
    
    double meancu;
    double var;
    
    /// ===== read the index file and create the molecules
    IndexFileReader   Index(indexfilename);
    std::vector<IndexGroup *> ndxgroups = Index.GetGroupIndex();
    
    int headid = 0;
    int wallid = 0;

    std::cout<<" We are reading the index file for molecules which flipflops \n";
    int i=0;
    std::cout<<"--------------------------------------------- "<<std::endl;
    std::cout<<"---The provided index file contains below groups ---- "<<std::endl;
    std::cout<<"--------------------------------------------- "<<std::endl;
    std::cout<<"--------------------------------------------- "<<std::endl;
    std::cout<<"-       group ID   |   group name - "<<std::endl;
    std::cout<<"--------------------------------------------- "<<std::endl;
    for (std::vector<IndexGroup *>::iterator it = ndxgroups.begin() ; it != ndxgroups.end(); ++it)
    {
        
        std::cout<<"        ("<<i<<")            "<<(*it)->GetGroupName()<<"  "<<std::endl;
        i++;
    }
    while (true)
    {
        std::cout<<" Please select a group for the molecules heads \n";
        std::cin>>headid;
        
        if(headid>=ndxgroups.size())
            std::cout<<" Error, the selected id is not accepted, please select again \n";
        else
            break;
    }
    while (true)
    {
        std::cout<<" Please select a group for the wall beads \n";
        std::cin>>wallid;
        
        if(wallid>=ndxgroups.size())
            std::cout<<" Error, the selected id is not accepted, please select again \n";
        else
            break;
    }
    IndexGroup * Group_Head = ndxgroups.at(headid);
    std::vector<int> pVHeadIndex = Group_Head->GetBeadsIndex();
    
    IndexGroup * Group_Wall = ndxgroups.at(wallid);
    std::vector<int> pVWallIndex = Group_Wall->GetBeadsIndex();
    bool systemhealth = true;

    if(pVWallIndex.size()!=m_WallPoints.size())
    {
        std::cout<<"Error: The head and tail group size is not equal, this means that you have selected wrong groups \n";
        systemhealth = false;
        
    }
    
    std::vector<bead*> pHeadBeads;

    for (int f=0;f<pVWallIndex.size();f++)
    {
        bead* b = m_pAllBeads.at(pVWallIndex.at(f)-1);
        m_pWallBeads.push_back(b);
        
    }
    for (int f=0;f<pVHeadIndex.size();f++)
    {
        bead* b = m_pAllBeads.at(pVHeadIndex.at(f)-1);
        pHeadBeads.push_back(b);
        
    }
    /////========= here we read a gro file and define a box and create the atoms
    

    ////
    int n=0;
    for (std::vector<bead* >::iterator it = m_pWallBeads.begin() ; it != m_pWallBeads.end(); ++it)
    {
        (*it)->UpdateC((m_WallPoints.at(n)).GetCurvature());
        n++;
    }
    
    ///============================================
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
// -------------------------------------------------- Time to do same with XTC file ---------------------------------------------------------
    if( trjfilename.at(trjfilename.size()-4)=='.')
    {
        if(trjfilename.at(trjfilename.size()-3)=='x' && trjfilename.at(trjfilename.size()-2)=='t' && trjfilename.at(trjfilename.size()-1)=='c' )
        {}
        else
        std::cout<<" the format of the trjectory file is not supported, please use xtc format "<<std::endl;

    }
    else
    {
        trjfilename=trjfilename+".xtc";
    }
    

                                      int *natoms;
                                      const char *file=trjfilename.c_str();
                                      const char *m="r";    /// the file should be read so mode is r
	                              XDRFILE *xdrfile1=xdrfile_open(file,  m);

if(xdrfile1 == NULL )
{
std::cout<<"\n"<<" are you kidding me, the xtc file does not exist: sorry we can not help you !!!"<<"\n";
}
else if(systemhealth == true)
{
//--------------------------------------------------------------------------------------------------------------------------
//------------------------------- if this part is exacuted means your xtc file exist so we do analysis ---------------------
//-----------------------------------------------------------------------------------------------------------------------------

    int count=0;
    int step(0),read(0);
    float time(1),prec(1000),lambda(1);
    matrix box;
    int natom;
    read_xtc_natoms(file,&natom);
    int oldstep=-1;
    bool splitwall= false;
    m_DX = 0.3;
 while (true)
 {// Loop 1
     

     rvec* xyz = new rvec[natom];
    
   read_xtc(xdrfile1,natom,&step,&time,box,&xyz[0],&prec);
     (*m_pBox)(0)=box[0][0];
     (*m_pBox)(1)=box[1][1];
     (*m_pBox)(2)=box[2][2];
     

     
     
     time=time/1000.0;

     
     //  for normal trajectory this should not be included
  //   if(step==0 && count>0)
   //      break;
     
     
     count++;
     
     if(step==oldstep)
         break;
     else
     oldstep=step;


     
     // Melanie does not like this
    // time=time-frame0time;
	 if (time>finaltime)
         break;

    // std::cout<<count<<" time "<<time<<"\n";

	if (time>=inialtime)
	NoFrame++;

     if(time >=inialtime && time<=finaltime && (NoFrame-1)%framejump==0)
     {
         acceptedframe++;

         if(acceptedframe%50==0)
          std::cout<<"time "<<time<<" step "<<step<<"\n";
         
             int j=0;
             for (std::vector<bead* >::iterator it = m_pAllBeads.begin() ; it != m_pAllBeads.end(); ++it)
             {
                 double x= xyz[j][0];
                 double y= xyz[j][1];
                 double z= xyz[j][2];
                 (*it)->UpdatePos(x, y, z);
                 j++;

             }
         if(splitwall==false)
         {
             FindUpperLowerWall();

             splitwall = true;
         }
         MakeSurfCurvatureHist();

         std::ofstream File1;
         File1.open("Upper.xvg");
         std::ofstream File3;
         File3.open("UpperF.xvg");
         std::ofstream File2;
         File2.open("Inner.xvg");
         
         for (i=0;i<1000;i++)
         {
             double pi=acos(-1);
             double h=3.7/2;
             double x0=2*pi/(*m_pBox)(0);
             double A=3;
             double t = double(i)*(*m_pBox)(0)/1000;
             double N=sqrt(1+A*A*x0*x0*cos(x0*t)*cos(x0*t));
             double x=t-h*A*x0*cos(x0*t)/N;
             double x2=t+h*A*x0*cos(x0*t)/N;

             if(x>(*m_pBox)(0))
                 x=x-(*m_pBox)(0);
             if(x<0)
                 x=x+(*m_pBox)(0);
             
             if(x2>(*m_pBox)(0))
                 x2=x2-(*m_pBox)(0);
             if(x2<0)
                 x2=x2+(*m_pBox)(0);
             
             
             double z0=A*sin(x0*t);
             double z=A*sin(x0*t)+ h/N; ;
             double z2=A*sin(x0*t) -h/N; ;

             File3<<t<<"  "<<x<<"   "<<z<<"   "<<x2<<"  "<<z2<<"  "<<z0<<"\n";
         }
         int s=0;
         for (std::vector<bead* >::iterator it = m_pWallBeads.begin() ; it != m_pWallBeads.end(); ++it)
         {
             Vec3D P = (*it)->GetPos();
             File1<<P(0)<<"   "<<P(2)<<"    "<<((m_WallPoints.at(s)).GetCurvature()).at(0)<<"\n";
             s++;
         }
         
         for (std::vector<bead* >::iterator it = m_pInWallBeads.begin() ; it != m_pInWallBeads.end(); ++it)
         {
             Vec3D P = (*it)->GetPos();
             File2<<P(0)<<"   "<<P(2)<<"\n";
         }

         /*
         double meancu = 0;
         double var = 0;
         double HeadSize = double(pHeadBeads.size());
         for (std::vector<bead* >::iterator it = pHeadBeads.begin() ; it != pHeadBeads.end(); ++it)
         {
             std::vector <double> C = FindCurv(*it);
             double c=C.at(0);
             
             meancu+= c/HeadSize;
             var+= c*c/HeadSize;
         }
         var = var-meancu*meancu;

         timeseri<<time<<"   "<<meancu<<"  "<<var<<"\n";
*/
         
         
        // GUN.Generate();



         

         

         
     }
     
     delete [] xyz;

 } // end Loop 1
    
    std::ofstream File4;
    File4.open("InnerF.xvg");
    for (int i=0;i<m_UpBeadCurve.size();i++)
    {
        
        File4<<(i+0.5)*m_DX<<"   "<<1/m_UpBeadCurve.at(i)<<"   "<<m_SUpBeadCurve.at(i)-(m_UpBeadCurve.at(i))*(m_UpBeadCurve.at(i))<<"   "<<m_LowerBeadCurve.at(i)<<"   "<<m_SLowerBeadCurve.at(i)-(m_LowerBeadCurve.at(i))*(m_LowerBeadCurve.at(i))<<"\n";
        
    }
}
    

    
}
Sin_1D::~Sin_1D()
{
    
}
std::vector<double> Sin_1D::FindCurv (bead *p)
{
    std::vector<double> C;
    double dist2 = 20000;
    int i=0;
    int idclosest = 0;
    for (std::vector<bead* >::iterator it = m_pWallBeads.begin() ; it != m_pWallBeads.end(); ++it)
    {
        
        double d2 = dist2between2bead((*it),p);
        if(d2<dist2)
        {
         dist2=d2;
         idclosest = i;
        }
        i++;
    }
    
    if(dist2>4)
        std::cout<<"warnning: best distance is larger than 4 nm \n";
    
    
    point po = m_WallPoints.at(idclosest);
    C = po.GetCurvature();
    return C;
    
}

double Sin_1D::dist2between2bead(bead* b1,bead* b2)
{
    
    double dist2=0;
    
    double x1=b1->GetXPos();
    double y1=b1->GetYPos();
    double z1=b1->GetZPos();
    
    double x2=b2->GetXPos();
    double y2=b2->GetYPos();
    double z2=b2->GetZPos();
    
    
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
void Sin_1D::FindUpperLowerWall()
{
    m_pUpWallBeads.clear();
    m_pInWallBeads.clear();
    int Nbins = int((*m_pBox)(0)/3);
    double binl =(*m_pBox)(0)/double(Nbins);
    
    
    for (int i=0;i<Nbins;i++)
    {
        std::vector<bead*> TemB;

        for (std::vector<bead* >::iterator it = m_pWallBeads.begin() ; it != m_pWallBeads.end(); ++it)
        {
            double X = (*it)->GetXPos();
            if(X>=i*binl && X<(i+1)*binl)
            {
                TemB.push_back(*it);
            }
        }
        double minZ = (*m_pBox)(2);
        double maxZ = 0;

        for (std::vector<bead* >::iterator it = TemB.begin() ; it != TemB.end(); ++it)
        {
            double Z = (*it)->GetZPos();
            if(minZ>Z)
            minZ = Z;
            if(maxZ<Z)
            maxZ = Z;
        }
        double midZ = (maxZ+minZ)/2;
        
        for (std::vector<bead* >::iterator it = TemB.begin() ; it != TemB.end(); ++it)
        {
            double Z = (*it)->GetZPos();
          //  if(Z>midZ)
            m_pUpWallBeads.push_back(*it);
          //  else if(Z<=midZ)
            m_pInWallBeads.push_back(*it);
         //   else
          //  std::cout<<"error1233: should not happen "<<midZ<<"   "<<Z<<"\n";
        }
        
    }
    
    
    std::cout<<"Note ----> among the wall beads, "<<m_pUpWallBeads.size()<<" is for upper monolayer and "<<m_pInWallBeads.size()<<" is for inner monolayer \n";
    
        std::cout<<"Note ----> all wall bead is "<<m_pWallBeads.size()<<"  and "<<m_pInWallBeads.size()+m_pUpWallBeads.size()<<" is sum of inn and upper \n";
    
    
    
    
    
}
void Sin_1D::MakeSurfCurvatureHist()
{


    int N=(*m_pBox)(0)/m_DX;
    m_DX = (*m_pBox)(0)/double(N);
    
    
    if(m_UpBeadCurve.size()==0)
    for (int i=0;i<N;i++)
    {
        m_UpBeadCurve.push_back(0);
        m_SUpBeadCurve.push_back(0);
        m_LowerBeadCurve.push_back(0);
        m_SLowerBeadCurve.push_back(0);
        m_NUp.push_back(0);
        m_NDown.push_back(0);
    }
    

for (int i=0;i<N;i++)
{
    for (std::vector<bead* >::iterator it = m_pUpWallBeads.begin() ; it != m_pUpWallBeads.end(); ++it)
    {
        Vec3D P = (*it)->GetPos();
        if(P(0)>=double(i)*m_DX && P(0)<double(i+1)*m_DX)
        {
            m_UpBeadCurve.at(i) = m_UpBeadCurve.at(i)+1/((*it)->GetCurvature()).at(0);
            m_SUpBeadCurve.at(i) = m_SUpBeadCurve.at(i)+(((*it)->GetCurvature()).at(0))*(((*it)->GetCurvature()).at(0));
            m_NUp.at(i) = m_NUp.at(i)+1;
        }
    }
}
    
    for (int i=0;i<N;i++)
    {
        for (std::vector<bead* >::iterator it = m_pInWallBeads.begin() ; it != m_pInWallBeads.end(); ++it)
        {
            Vec3D P = (*it)->GetPos();
            if(P(0)>=double(i)*m_DX && P(0)<double(i+1)*m_DX)
            {
                m_LowerBeadCurve.at(i) = m_LowerBeadCurve.at(i)+1/((*it)->GetCurvature()).at(0);
                m_SLowerBeadCurve.at(i) = m_SLowerBeadCurve.at(i)+(((*it)->GetCurvature()).at(0))*(((*it)->GetCurvature()).at(0));
                m_NDown.at(i) = m_NDown.at(i)+1;
            }
        }
    }
    
    for (int i=0;i<N;i++)
    {
 
        m_LowerBeadCurve.at(i) = m_LowerBeadCurve.at(i)/m_NDown.at(i);
        m_SLowerBeadCurve.at(i) = m_SLowerBeadCurve.at(i)/m_NDown.at(i);

        m_UpBeadCurve.at(i) = m_UpBeadCurve.at(i)/m_NUp.at(i);
        m_SUpBeadCurve.at(i) = m_SUpBeadCurve.at(i)/m_NUp.at(i);

  
    }
}




#endif



