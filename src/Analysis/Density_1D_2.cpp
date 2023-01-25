 #if !defined(AFX_Density_1D_2_CPP_7F4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_Density_1D_2_CPP_7F4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "Density_1D_2.h"
#include "GroFile.h"
#include "IndexFileReader.h"
#include  "IndexGroup.h"
#include  "GenerateUnitCells.h"
#include "PDBFile.h"


Density_1D_2::Density_1D_2(Argument *pArgu)
{


//////////  Atom number in the input file is the number of the atom in the coordinate file and it start from 1.... but index which is used in reading xtc started from 0
///////// because that we use Atomno-1

    //std::string pdbfilename = pArgu->GetpdbFileName();
    //PDBFile pdb;
    //m_WallPoints = pdb.ReadPDBFile(pdbfilename);
    // GenerateUnitCells   GUN(m_pAllBeads,pArgu,m_pBox);

std::cout<<"\n";
std::cout<<"------------- Density_1D_2 is running: Which mean your inputs were ok, but let see if its make sense too  -------------"<<"\n";
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

    LMatrix RefMat(Nx,2);
    LMatrix TargetMat(Nx,2);
    LMatrix W(Nx,2);

    

    
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
         LMatrix M = SplitDenisty(TemTargetMat,TemRefMat);

         for (int i=0;i<Nx;i++)
         {
             for (int j=0;j<Nz;j++)
             {

                 if(M(i,3)!=0)
                 {
                     TargetMat(i,0) = TargetMat(i,0)+M(i,1)/M(i,3);
                     W(i,0)=W(i,0)+1;
                 }
                 if(M(i,4)!=0)
                 {
                     TargetMat(i,1) = TargetMat(i,1)+M(i,2)/M(i,4);
                     W(i,1)=W(i,1)+1;
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

            if(W(i,0)!=0)
            {
                
                TargetMat(i,0) = TargetMat(i,0)/W(i,0)/double(pTargetBeads.size())*double(pReferenceBeads.size());
            }
        if(W(i,1)!=0)
        {
            
            TargetMat(i,1) = TargetMat(i,1)/W(i,1)/double(pTargetBeads.size())*double(pReferenceBeads.size());
        }

    }


    ///******* end rescale with respect to accepted frames

    

}
    
    std::string gfilename = pArgu->GetGeneralOutputFilename();
    std::string outputfilename=gfilename+targetname+"-"+referencename+".xvg";
    std::ofstream output;
    output.open(outputfilename.c_str());
    
    output<<"#####   i     x      "<<targetname<<"updensity   "<<targetname<<"downdensity   "<<referencename<<"updensity   "<<referencename<<"downdensity   midsurfcetarget midsurfaceref "<<"\n";
    for (int i=0;i<Nx;i++)
    {
        output<<i<<"   "<<BLx*i<<"   "<<TargetMat(i,0)<<"   "<<TargetMat(i,1)<<"\n";
    }
    

    
}
Density_1D_2::~Density_1D_2()
{
    
}
LMatrix Density_1D_2::SplitDenisty(LMatrix a, LMatrix ref)
{
    int Ny = a.Csize();
    int Nx = a.Rsize();
    
    LMatrix SplitMat(Nx,5);

    for (int i=0;i<Nx;i++)
    {
        int ymax = 0;
        int ymin = 0;
        double upper = 0;
        double inner = 0;
        double uref = 0;
        double iref = 0;
        for (int j=0;j<Ny;j++)
        {
            if(ref(i,j)!=0 && ymin==0)
                ymin = j;
            else if(ref(i,j)!=0)
                ymax = j;
        }
        int mid = int(double(ymax+ymin)/2.0);
        for (int j=0;j<Ny;j++)
        {
            if(j>mid)
            {
                upper+= a(i,j);
                uref+=ref(i,j);
            }
            else
            {
                inner+= a(i,j);
                iref+=ref(i,j);
            }

        }
        if(double(ymax-ymin)/double(Ny)*16.95757>8 ||double(ymax-ymin)/double(Ny)*16.95757<2)
        {
        std::cout<<"  not expected "<<ymax-ymin<<"   "<<Ny<<"   "<<double(ymax-ymin)/double(Ny)*16.95757<<"\n";
            SplitMat(i,0) = 0;
            SplitMat(i,1) = 0;
            SplitMat(i,2) = 0;
            SplitMat(i,3) = 0;
            SplitMat(i,4) = 0;
        }
        else
        {
        SplitMat(i,0) = mid;
        SplitMat(i,1) = upper;
        SplitMat(i,2) = inner;
        SplitMat(i,3) = uref;
        SplitMat(i,4) = iref;
        }
    }
    
    
    return SplitMat;

}




#endif



