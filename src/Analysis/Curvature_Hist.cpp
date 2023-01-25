 #if !defined(AFX_Curvature_Hist_CPP_7F4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_Curvature_Hist_CPP_7F4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "Curvature_Hist.h"
#include "GroFile.h"
#include "IndexFileReader.h"
#include  "IndexGroup.h"
#include  "GenerateUnitCells.h"
#include "PDBFile.h"

Curvature_Hist::Curvature_Hist(Argument *pArgu)
{


//////////  Atom number in the input file is the number of the atom in the coordinate file and it start from 1.... but index which is used in reading xtc started from 0
///////// because that we use Atomno-1



std::cout<<"\n";
std::cout<<"------------- Curvature_Hist is running: Which mean your inputs were ok, but let see if its make sense too  -------------"<<"\n";
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
    
    timeseri<<"@@     time (ns)    <H>  <H*H>     <K>    <K*K>"<<"\n";
    
    
    double meancu;
    double var;
    
    
    //====== Wall Histogram
   // m_WallPoints
    std::vector<double> W_C;
    std::vector<double> W_K;
    std::vector<double> L_C;
    std::vector<double> L_K;
    double hmin = 10000000;
    double hmax = -10000000;
    double kmin = 10000000;
    double kmax = -10000000;
    for (std::vector<point >::iterator it = m_WallPoints.begin() ; it != m_WallPoints.end(); ++it)
    {
       std::vector<double> C = (*it).GetCurvature();
        W_C.push_back(C.at(0));
        W_K.push_back(C.at(1));
        
        if(C.at(0)>hmax)
        hmax = C.at(0);
        else if(C.at(0)<hmin)
        hmin = C.at(0);
        
        if(C.at(1)>kmax)
        kmax = C.at(1);
        else if(C.at(1)<kmin)
        kmin = C.at(1);
    }


    

    
    


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
         double meancu = 0;
         double var = 0;
         double gauss = 0;
         double Gvar = 0;
         double HeadSize = double(pHeadBeads.size());
         for (std::vector<bead* >::iterator it = pHeadBeads.begin() ; it != pHeadBeads.end(); ++it)
         {
             std::vector <double> C = FindCurv(*it);
             double  c=C.at(0);
             double  K=C.at(1);
             L_C.push_back(c);
             L_K.push_back(K);
             meancu+= c/HeadSize;
             var+= c*c/HeadSize;
             gauss+=K/HeadSize;
             Gvar+=K*K/HeadSize;
             
             
             
         }
         var = var-meancu*meancu;
         
         
         
         

         timeseri<<time<<"   "<<meancu<<"  "<<var<<"   "<<gauss<<"   "<<Gvar<<"\n";

         
         
        // GUN.Generate();

         if(isfirstframe == true)
         {
             isfirstframe =false;
             frame0time = time;

         }

         

         

         
     }
     
     delete [] xyz;

 } // end Loop 1
}
    
    
    std::ofstream OS_C_MAT;
    OS_C_MAT.open("C_MAT.xvg");
    double dx = 0.001;
    int NH = (hmax -hmin)/dx;
    int NK = (kmax -kmin)/dx;
    
    

    double C_MAT[NH][NK];
    double WG = double(L_C.size());
    for (int i = 0 ; i<L_C.size();i++)
    {
        for (int n = 0 ; n<NH;n++)
        for (int m = 0 ; m<NK;m++)
        {
            if (double(n)*dx+hmin <L_C.at(i) && double(n+1)*dx+hmin >L_C.at(i))
            if (double(m)*dx+kmin <L_K.at(i) && double(m+1)*dx+kmin >L_K.at(i))
            C_MAT[n][m] = C_MAT[n][m]+1.0/WG;
        }
        
    }
    
    for (int n = 0 ; n<NH;n++)
    {
        for (int m = 0 ; m<NK;m++)
        {
        
            OS_C_MAT<<C_MAT[n][m] <<"     ";

        }
        OS_C_MAT<<"\n";
    }
    
    std::vector<std::vector<double> > W_H_C = hist(W_C);
    std::vector<std::vector<double> > W_H_K = hist(W_K);
    std::vector<std::vector<double> > L_H_C = hist(L_C);
    std::vector<std::vector<double> > L_H_K = hist(L_K);
    
    
    std::ofstream hist_H;
    hist_H.open("OS_H_HIST.xvg");
    
    std::ofstream hist_K;
    hist_K.open("OS_K_HIST.xvg");
    
    int Size = (W_H_C.at(0)).size();
    for (int i = 0; i<Size;i++)
    {
        if((L_H_C.at(0)).size()>i)
           hist_H<<(W_H_C.at(0)).at(i)<<"    "<<(W_H_C.at(1)).at(i)<<"    "<<(L_H_C.at(0)).at(i)<<"    "<<(L_H_C.at(1)).at(i)<<"\n";
        else
           hist_H<<(W_H_C.at(0)).at(i)<<"    "<<(W_H_C.at(1)).at(i)<<"    -          -     "<<"\n";

    }
           
    Size = (W_H_K.at(0)).size();
    for (int i = 0; i<Size;i++)
    {
        if((L_H_C.at(0)).size()>i)
            hist_K<<(W_H_K.at(0)).at(i)<<"    "<<(W_H_K.at(1)).at(i)<<"    "<<(L_H_K.at(0)).at(i)<<"    "<<(L_H_K.at(1)).at(i)<<"\n";
        else
            hist_K<<(W_H_K.at(0)).at(i)<<"    "<<(W_H_K.at(1)).at(i)<<"    -          -     "<<"\n";
                  
    }
    
    
    

    
}
Curvature_Hist::~Curvature_Hist()
{
    
}
std::vector<double> Curvature_Hist::FindCurv (bead *p)
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

double Curvature_Hist::dist2between2bead(bead* b1,bead* b2)
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
std::vector<std::vector<double> > Curvature_Hist::hist(std::vector<double> data)
{
    std::vector<double> p;
    std::vector<double> x;
    std::vector<std::vector<double> > h;
    double dl = 0.001;
    double xmin = 10000000;
    double xmax = -10000000;
    for (int i=0;i<data.size();i++)
    {
        if(data.at(i)>xmax)
        xmax = data.at(i);
        else if(data.at(i)<xmin)
        xmin = data.at(i);
    }
    int bin = (xmax - xmin)/dl+1;
    double Hist[bin][2];

    for (int j=0;j<bin;j++)
    {
        Hist[j][1] = 0;
        Hist[j][0] = 0;
        
    }
    
    double x0= (xmax - xmin);
    double dx = (xmax - xmin)/bin;
    double tot = double (data.size());
    
    for (int i=0;i<data.size();i++)
    {
        for (int j=0;j<bin;j++)
        {
            if(data.at(i)<xmin+dx*double(j+1) && data.at(i)>xmin+dx*double(j))
            {
                Hist[j][0] = dx*(double(j)+0.5);
                Hist[j][1] = Hist[j][1] + 1.0/tot;
                
            }
        }
    }
    
    for (int j=0;j<bin;j++)
    {
        p.push_back(Hist[j][1]);
        x.push_back(Hist[j][0]);
        
    }
    
    h.push_back(x);
    h.push_back(p);
    
    return h;
}



#endif



