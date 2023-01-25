 #if !defined(AFX_Density_1D_CPP_7F4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_Density_1D_CPP_7F4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "Density_1D.h"
#include "GroFile.h"
#include "IndexFileReader.h"
#include  "IndexGroup.h"
#include  "GenerateUnitCells.h"
#include "PDBFile.h"


Density_1D::Density_1D(Argument *pArgu)
{


//////////  Atom number in the input file is the number of the atom in the coordinate file and it start from 1.... but index which is used in reading xtc started from 0
///////// because that we use Atomno-1

    //std::string pdbfilename = pArgu->GetpdbFileName();
    //PDBFile pdb;
    //m_WallPoints = pdb.ReadPDBFile(pdbfilename);
    // GenerateUnitCells   GUN(m_pAllBeads,pArgu,m_pBox);

std::cout<<"\n";
std::cout<<"------------- Density_1D is running: Which mean your inputs were ok, but let see if its make sense too  -------------"<<"\n";
    std::cout<<"------------- Please note: time unit is ns in this software  -------------"<<"\n";

std::cout<<"\n";


Nfunction f;      // In this class there are some useful function and we can use it.


    
    
    
    //***************************************
    //********* read the index file
    std::string indexfilename = pArgu->GetIndexFileName();
    IndexFileReader   Index(indexfilename);
    std::vector<IndexGroup *> ndxgroups = Index.GetGroupIndex();
    std::string targetname  = pArgu->GetTargetGroupName();
    std::string referencename = pArgu->GetReferenceGroupName();
    bool tg = false;
    bool rg = false;
    int refGid = 0;
    int targetGid = 0;
    int GroupID = 0;
    bool m_health = true;
    for (std::vector<IndexGroup *>::iterator it = ndxgroups.begin() ; it != ndxgroups.end(); ++it)
    {
        
        if((*it)->GetGroupName()==targetname)
        {
            tg = true;
            targetGid = GroupID;

            
        }
        if((*it)->GetGroupName()==referencename)
        {
            rg = true;
            refGid = GroupID;
        }
        GroupID++;
    }
    if(rg!=true  || tg!=true)
    {
        m_health = false;
        std::cout<<"error--->ref or traget groups does not exist in the index file "<<std::endl;
        std::cout<<"--------------------------------------------- "<<std::endl;
        std::cout<<"---The provided index file contains below groups ---- "<<std::endl;
        std::cout<<"--------------------------------------------- "<<std::endl;
        std::cout<<"--------------------------------------------- "<<std::endl;
        std::cout<<"-       group ID   |   group name - "<<std::endl;
        std::cout<<"--------------------------------------------- "<<std::endl;
        int i =0;
        for (std::vector<IndexGroup *>::iterator it = ndxgroups.begin() ; it != ndxgroups.end(); ++it)
        {
        
            std::cout<<"        ("<<i<<")            "<<(*it)->GetGroupName()<<"  "<<std::endl;
            i++;
        }
    }
    
    IndexGroup * Group_Target = ndxgroups.at(targetGid);
    std::vector<int> pVTargetIndex = Group_Target->GetBeadsIndex();
    IndexGroup * Group_Reference = ndxgroups.at(refGid);
    std::vector<int> pVReferenceIndex = Group_Reference->GetBeadsIndex();
    //***************************************
    //********* end read the index file
    
    //***************************************
    //********* gro file
    std::string grofilename = pArgu->GetGroFileName();
    GroFile Gro(grofilename);
    m_pBox = Gro.GetBox();
    m_pAllBeads = Gro.GetpAllBeads();
    std::cout<<"  System contains "<<m_pAllBeads.size()<<"   beads \n";
    std::cout<<"  box info "<<(*m_pBox)(0)<<"  "<<(*m_pBox)(1)<<"  "<<(*m_pBox)(2)<<"   beads \n";

    //************************************************************
    //***** find the beads of target and reference beads
    //**********************************************************************
    std::vector<bead*> pTargetBeads;
    std::vector<bead*> pReferenceBeads;

    for (int f=0;f<pVReferenceIndex.size();f++)
    {
        bead* b = m_pAllBeads.at(pVReferenceIndex.at(f)-1);
        pReferenceBeads.push_back(b);
        
    }
    for (int f=0;f<pVTargetIndex.size();f++)
    {
        bead* b = m_pAllBeads.at(pVTargetIndex.at(f)-1);
        pTargetBeads.push_back(b);
        
    }
    std::cout<<" target group has "<<pTargetBeads.size()<<" beads and ref group has "<<pReferenceBeads.size()<<" beads \n";
    //******** ref mat and target mat ************************************************************************
    double BL = pArgu->GetBinLength();
    int Nx = int((*m_pBox)(0)/BL);
    int Nz = int((*m_pBox)(2)/BL);
    double BLx =(*m_pBox)(0)/double(Nx);
    double BLz =(*m_pBox)(2)/double(Nz);

    LMatrix RefMat(Nx,Nz);
    LMatrix TargetMat(Nx,Nz);
    LMatrix W(Nx,Nz);

    

    
    //******** end ref mat and target mat ************************************************************************

    //********************************************************************************
    //******** Trajectory file name ************************************************************************
    int framejump = pArgu->GetJump();
    int NoFrame = 0;
    int acceptedframe = 0;
    std::string trjfilename = pArgu->GetTrjFile();
    /// read data from argument class
    int inialtime = pArgu->GetBStep();
    int finaltime = pArgu->GetEStep();
    
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
else if(m_health == true)
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
   /* std::ofstream out;
    out.open("test.xvg");*/
 while (true)
 {// Loop 1
     

     rvec* xyz = new rvec[natom];
    
   read_xtc(xdrfile1,natom,&step,&time,box,&xyz[0],&prec);
     (*m_pBox)(0)=box[0][0];
     (*m_pBox)(1)=box[1][1];
     (*m_pBox)(2)=box[2][2];
     

    // std::cout<<"  box info "<<(*m_pBox)(0)<<"  "<<(*m_pBox)(1)<<"  "<<(*m_pBox)(2)<<"   beads \n";

     
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
         
             int beadidc=0;
             for (std::vector<bead* >::iterator it = m_pAllBeads.begin() ; it != m_pAllBeads.end(); ++it)
             {
                 double x= xyz[beadidc][0];
                 double y= xyz[beadidc][1];
                 double z= xyz[beadidc][2];
                 (*it)->UpdatePos(x, y, z);
                 beadidc++;

             }

         LMatrix TemRefMat(Nx,Nz);
         LMatrix TemTargetMat(Nx,Nz);

            for (std::vector<bead *>::iterator it = pTargetBeads.begin() ; it != pTargetBeads.end(); ++it)
            {
                double X = (*it)->GetXPos();
                double Z = (*it)->GetZPos();
                //out<<X<<"  "<<Y<<"\n";
               // std::cout<<(*it)->GetResName()<<"  "<<(*it)->GetBeadName()<<"\n";
                 
                 int i = int (X*double(Nx)/(*m_pBox)(0));
                 int j = int (Z*double(Nz)/(*m_pBox)(2));
                if(i>  Nx-1 || j>Nz-1)
                {
                    std::cout<<X<<"   "<<Nx<<"  "<<(*m_pBox)(0)<<"\n";
                    std::cout<<Z<<"   "<<Nz<<"  "<<(*m_pBox)(2)<<"\n";

                }
                TemTargetMat(i,j) = TemTargetMat(i,j)+1;
                // TargetMat(i,j)= TargetMat(i,j)+1.0/double(pTargetBeads.size());
                
            }
             for (std::vector<bead *>::iterator it = pReferenceBeads.begin() ; it != pReferenceBeads.end(); ++it)
             {
                 double X = (*it)->GetXPos();
                 double Z = (*it)->GetZPos();
                 //out<<X<<"  "<<Y<<"\n";
                 // std::cout<<(*it)->GetResName()<<"  "<<(*it)->GetBeadName()<<"\n";
                 
                 int i = int (X*double(Nx)/(*m_pBox)(0));
                 int j = int (Z*double(Nz)/(*m_pBox)(2));
                 if(i>  Nx-1 || j>Nz-1)
                 {
                     std::cout<<X<<"   "<<Nx<<"  "<<(*m_pBox)(0)<<"\n";
                     std::cout<<Z<<"   "<<Nz<<"  "<<(*m_pBox)(2)<<"\n";
                     
                 }
                  TemRefMat(i,j)= TemRefMat(i,j)+1.0;
             }
         for (int i=0;i<Nx;i++)
         {
             for (int j=0;j<Nz;j++)
             {
                 if(TemTargetMat(i,j)!=0)
                 {
                     TargetMat(i,j) = TargetMat(i,j)+TemTargetMat(i,j)/TemRefMat(i,j);
                     W(i,j)=W(i,j)+1;
                 }
                 else if(TemRefMat(i,j)!=0)
                 {
                     W(i,j)=W(i,j)+1;
                 }
                 // RefMat(i,j) = RefMat(i,j)/double(acceptedframe);
                 
             }
         }
     }

     
     delete [] xyz;

 } // end Loop 1
    ///******* rescale with respect to accepted frames
    for (int i=0;i<Nx;i++)
    {
        for (int j=0;j<Nz;j++)
        {
            if(TargetMat(i,j)!=0)
            {TargetMat(i,j) = TargetMat(i,j)/W(i,j)/double(pTargetBeads.size())*double(pReferenceBeads.size());
            RefMat(i,j) = RefMat(i,j)/W(i,j);
            }

        }
    }
    {
        std::string outputfilename=targetname+"mat.dat";
        std::ofstream output;
        output.open(outputfilename.c_str());
        for (int i=0;i<Nx;i++)
        {
            for (int j=0;j<Nz;j++)
            {
                output<<TargetMat(i,j)<<"   ";
                
            }
            output<<"\n";
        }
    }
    {
        std::string outputfilename=referencename+"mat.dat";
        std::ofstream output;
        output.open(outputfilename.c_str());
        for (int i=0;i<Nx;i++)
        {
            for (int j=0;j<Nz;j++)
            {
                output<<RefMat(i,j)<<"   ";
                
            }
            output<<"\n";
        }
    }
    ///******* end rescale with respect to accepted frames

    

}
    
    std::string gfilename = pArgu->GetGeneralOutputFilename();
    LMatrix RefVec= SplitDenisty(RefMat);
    LMatrix TagetVec= SplitDenisty(TargetMat);
    std::string outputfilename=gfilename+targetname+"-"+referencename+".xvg";
    std::ofstream output;
    output.open(outputfilename.c_str());
    
    output<<"#####   i     x      "<<targetname<<"updensity   "<<targetname<<"downdensity   "<<referencename<<"updensity   "<<referencename<<"downdensity   midsurfcetarget midsurfaceref "<<"\n";
    for (int i=0;i<Nx;i++)
    {
        output<<i<<"   "<<BLx*i<<"   "<<TagetVec(i,1)<<"   "<<TagetVec(i,2)<<"    "<<RefVec(i,1)<<"   "<<RefVec(i,2)<<"  "<<BLz*TagetVec(i,0)<<"  "<<BLz*RefVec(i,0)<<"\n";
    }
    

    
}
Density_1D::~Density_1D()
{
    
}
LMatrix Density_1D::SplitDenisty(LMatrix a)
{
    int Ny = a.Csize();
    int Nx = a.Rsize();
    
    std::cout<<Nx<<"  "<<Ny<<"\n";
    LMatrix SplitMat(Nx,3);

    for (int i=0;i<Nx;i++)
    {
        int ymax = 0;
        int ymin = 0;
        double upper = 0;
        double inner = 0;
        for (int j=0;j<Ny;j++)
        {
            if(a(i,j)!=0 && ymin==0)
                ymin = j;
            else if(a(i,j)!=0)
                ymax = j;
        }
        int mid = int(double(ymax+ymin)/2.0);
        for (int j=0;j<Ny;j++)
        {
            if(j>mid)
                upper+= a(i,j);
            else
                inner+= a(i,j);
        }
        if(double(ymax-ymin)/double(Ny)*16.95757>7 ||double(ymax-ymin)/double(Ny)*16.95757<3)
        std::cout<<"  not expected "<<ymax-ymin<<"   "<<Ny<<"   "<<double(ymax-ymin)/double(Ny)*16.95757<<"\n";
        SplitMat(i,0) = mid;
        SplitMat(i,1) = upper;
        SplitMat(i,2) = inner;
    }
    
    
    return SplitMat;

}




#endif



