  #if !defined(AFX_BackMap_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_BackMap_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "BackMap.h"
#include "GroFile.h"
#include "GenerateUnitCells.h"
#include "Def.h"
#include "PDBFile.h"
#include "GenDomains.h"
#include "PointBasedBlueprint.h"
/*
 1) read the point// checked!
 2) exclude the point base of exclsuion data// checked
 3) place proteins or generate them // checked
 4) exclude point based on proteins  // checked
 5) place the lipids
 6) check for all unnessry variable and function and remove them
 6) make protein placment better
 
 
 
 further extentions
 
 1). in random inc, we should prevent them to be close
 2.) angle and theta is not used yet
 3.) pattern based protein insertion
 
 
 
 */
BackMap::BackMap(Argument *pArgu)
{
    m_monolayer = false;

    srand (pArgu->GetSeed());
    std::cout<<"\n";
    std::cout<<"=========================================================================================================="<<"\n";
    std::cout<<"****************** "<< SoftWareName <<" ******************   "<<"\n";
    std::cout<<"******************  Version:  "<<SoftWareVersion<<" ****************** \n";
    std::cout<<"=========================================================================================================="<<"\n";
    Nfunction f;      // In this class there are some useful function and we can use it.

    //====== getting data points to create cg membrane
    std::cout<<"---> attempting to obtain point data \n";
    PointBasedBlueprint SurfDataPoint(pArgu);
    std::vector<point*>  pPointUp = SurfDataPoint.m_pPointUp;
    std::vector<point*>  pPointDown = SurfDataPoint.m_pPointDown;
    std::vector<exclusion*>  m_pExc = SurfDataPoint.m_pExc;
    std::vector<inclusion*>  m_pInc = SurfDataPoint.m_pInc;
    Vec3D *m_pBox = SurfDataPoint.m_pBox;
    Wall *m_pWall = SurfDataPoint.m_pWall;
    m_monolayer = SurfDataPoint.m_monolayer;
    std::cout<<"---> point data has been obtained \n";

    //======== OutPut file name declaration and finding input file names ========================================
    std::string gname = pArgu->GetGeneralOutputFilename();   // get the generic name for the outputs
    m_FinalOutputGroFileName =gname+".gro";                     // create the output gro file name based on the generic name
    std::string m_FinalTopologyFileName=gname+".top";
    
    //======
    std::cout<<"---> attempting to generate molecule type \n";
    GenerateMolType  MOLTYPE(pArgu);   // using the str file, the included gro file in the str and the lib file, different mol types will be generated. proteins and lipids are treated as mols.
    m_map_MolName2MoleculesType = MOLTYPE.GetMolType();   // a map containing all the mol types: (name, MolType)
    std::cout<<"---> molecule types have been generated \n";
  
    //== we should exclude points and get rid of exclusion. This is done by making the area of the point zero.
    //== this could be made more efficient but is not needed as exclusion should not inlcude to many points
    ExcludePointsUsingExclusion(m_pExc, pPointUp, pPointDown);
    
    //==== now we need to place the proteins
    // first we read str file to find protein info
    std::string strfilename = pArgu->GetStructureFileName();    // str file name
    if(FindProteinList(strfilename)==false) // this data will be stored in m_map_IncID2ProteinLists map(tsi_protein_id,ProteinList)
        std::exit(0);
    //== before going further, it would be better to check if the info about
    //proteins in str file also exist in the PLM output and also in the mol type
    if(CheckProteinInfo (m_map_IncID2ProteinLists, m_map_MolName2MoleculesType, m_pInc)==false) //
        std::exit(0);
    

    if(m_pInc.size()!=0)
    {
        //== I do not know why do we need this
        for ( std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin(); it != m_map_IncID2ProteinLists.end(); it++ )
        {
            int number = 0;
            int id = (it->second).ID;
            for ( std::vector<inclusion*>::iterator it1 = m_pInc.begin(); it1 != m_pInc.end(); it1++ )
            {
                if(id==(*it1)->GetTypeID())
                    number++;
            }
            (it->second).created = number;
        }
    }
    else  // if the inclusion file is empty, we generate proteins according to the data in the str file
    {
        CreateRandomInclusion(pPointUp);
        for (std::vector<inclusion>::iterator it2 = m_RandomInc.begin() ; it2 != m_RandomInc.end(); ++it2)
            m_pInc.push_back(&(*it2));
    }
    
//== better to place the proteins and make then remove the excluded points  then generate the domains. Note, we should tell the user some domain could be removed due to protein and exclsuion
    //== however, if a error happen due to user, the error will be revealed later about domain. We can add a check function later
    
    m_InclusionDirectionType = pArgu->GetInclusionDirectionType(); //Note: in the normal condition PLM always write global so only applicable if you want to change the point folder manually
    if(m_InclusionDirectionType=="Local")
    {
        std::cout<<"The inclusion direction type is: "<<m_InclusionDirectionType<<"\n";
        std::cout<<" Note: in the normal condition PLM always write global so only applicable if you want to change the point folder manually\n";
    }
    //== Placing the inclusions;
    if(PlaceProteins(pPointUp)==false) // this creates all the protein beads and put them in m_FinalBeads
        std::exit(0);
    std::cout<<"---> proteins are placed, now we remove points that are close to the proteins \n";
    
    //=== removing points closeby the proteins
    {
        double RCutOff = pArgu->GetRCutOff();     /// That will be counted as a cutoff for the protein_lipid distance
        std::vector<bead*> tempropbeads;
    
        for ( std::vector<bead>::iterator it = m_FinalBeads.begin(); it != m_FinalBeads.end(); it++ )
            tempropbeads.push_back(&(*it));
        
        bool RMpoint = RemovePointsCloseToBeadList(pPointUp, pPointDown, tempropbeads, RCutOff, m_pBox);
        if(RMpoint==false)
            std::exit(0);
    }
    std::cout<<"---> generating domains using the input files \n";


    //== next task is to make the domain work and also make sure that the point with area zero is not part of the domain points

    
    //std::string strfilename = pArgu->GetStructureFileName();    // Get the stracture file name
    //m_ResID = 1;
   //m_Renormalizedlipidratio = pArgu->GetRenorm();
   //m_Iter = pArgu->GetIter();
    // Perhaps reading the str file to find out which lipids are needed and which proteins are there ....
    //GenDomains GENDOMAIN(strfilename,pPointUp,pPointUp,m_Renormalizedlipidratio);  // this somehow reads the lipids
    //std::vector<Domain*> pAllDomain = GENDOMAIN.GetDomains();


    

    /*


//==========================================================================================================



    
    //============================== Finding the total area of the layers
    m_TotalAreaUp = 0.0;
    m_TotalAreaDown = 0.0;
    for ( std::vector<point*>::iterator it = m_point1.begin(); it != m_point1.end(); it++ )
        m_TotalAreaUp+= (*it)->GetArea();
    for ( std::vector<point*>::iterator it = m_point2.begin(); it != m_point2.end(); it++ )
        m_TotalAreaDown+= (*it)->GetArea();

    if(m_monolayer == false)
    std::cout<<"--> Note: the total upper monolayer area is "<<m_TotalAreaUp<<" and the total lower monolayer area is "<<m_TotalAreaDown<<"\n";
    else
    std::cout<<"--> Note: the total monolayer area is "<<m_TotalAreaUp<<" \n";


    }
    // Make all the domain containing different lipids
    GenDomains GENDOMAIN(strfilename,p1,p2,m_Renormalizedlipidratio);
    std::vector<Domain*> pAllDomain = GENDOMAIN.GetDomains();
    std::cout<<" Number of the domains defined in the input file  "<< pAllDomain.size()/2 <<"\n";
    int layer = 0;
        std::cout<<"------  we aim to generate  ------- \n";
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
    

//    for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
    {
        
    }

    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {
            std::vector<DomainLipid*> DL = (*it)->GetpDomainLipids();

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
                if (m_map_MolName2MoleculesType.count(ltype) == 0)
                    std::cout << "Error:-----> molecule name " <<ltype<<" does not exist in the lib files \n";
                
                GenLipid(m_map_MolName2MoleculesType.at(ltype), 0, Pos, N, Dir, T1, T2);
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
     */
    std::cout<<"*************************** Number of created Lipids,   ********************** \n";
    /*layer = 0;
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
    
    for ( std::map<int,ProteinList>::iterator it = m_TotalProteinList.begin(); it != m_TotalProteinList.end(); it++ )
        Topgro<<(it->second).ProteinName<<"   "<<(it->second).created<<"\n";

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

*/

}
BackMap::~BackMap()
{
    
}
void BackMap::GenLipid(MolType moltype, int listid, Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2)
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
void BackMap::GenProtein(MolType moltype, int listid, Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2)
{

        Tensor2 LG = TransferMatLG(Normal, t1, t2);
        Tensor2 GL = LG.Transpose(LG);


        //===== to fit to the protein diretion
        Vec3D LocalDir;
        if(m_InclusionDirectionType=="Global") //Note: in the normal condition PLM always write global so only applicable if you want to change the point folder manually
        LocalDir = GL*Dir;
        else if(m_InclusionDirectionType=="Local")
        LocalDir = Dir;
        double C=LocalDir(0);
        double S= LocalDir(1);
        Tensor2 Rot=Rz(C,S);
        Vec3D DH= Normal*((m_map_IncID2ProteinLists.at(listid)).Z0);
        double phi = (m_map_IncID2ProteinLists.at(listid)).Phi;
        double theta = (m_map_IncID2ProteinLists.at(listid)).Theta;

        //
        std::vector<bead> vbeads = moltype.Beads;
        for ( std::vector<bead>::iterator it = vbeads.begin(); it != vbeads.end(); it++ )
        {
            Vec3D BPos((*it).GetXPos(),(*it).GetYPos(),(*it).GetZPos());
            Vec3D vX = LG*(Rot*BPos)+ Pos+DH;
            int beadid = (m_FinalBeads.size()+1);
            bead TemB(beadid, (*it).GetBeadName(), (*it).GetBeadType(), (*it).GetResName(), m_ResID, vX(0), vX(1),vX(2));
            m_FinalBeads.push_back(TemB);
        }
        m_ResID++;
}
double BackMap::dist2between2Points(Vec3D X1,Vec3D X2)
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
void BackMap::WriteFinalGroFile()
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
        std::string rname  = (*it).GetResName();
        std::string bname = (*it).GetBeadName();
        
        const char* A1=rname.c_str();
        const char* A2=bname.c_str();
        int resid=(*it).GetResid();
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid%100000,A1,A2,i%100000,x,y,z );
        
    }
    
    
    fprintf(fgro,  "%10.5f%10.5f%10.5f\n",m_Box(0),m_Box(1),m_Box(2) );
    fclose(fgro);
    
    
}

Tensor2 BackMap::Rz(double cos, double sin)
{
    Tensor2  R;
    
    R(0,0) = cos;
    R(0,1) = -sin;
    R(1,0) = sin;
    R(1,1) = cos;
    R(2,2) = 1;

    
    
    return R;
}
//=== we can make this function better
void BackMap::CreateRandomInclusion(std::vector<point*> &pPointUp)
{
    int totinc = 0;
    int totcreated = 0;
    //==== reading input file to find number of requested proteins
    if(m_map_IncID2ProteinLists.size()!=0)
    {
        std::cout<<" According to the data and area we generate  "<<m_map_IncID2ProteinLists.size()<<" proteins types \n";
    }
    for ( std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin(); it != m_map_IncID2ProteinLists.end(); it++ )
    {
        (it->second).created = 0;
        double ratio = (it->second).Ratio;
        std::string type = (it->second).ProteinName;
        double area = (m_map_MolName2MoleculesType.at(type)).molarea;
        int neededno = m_TotalAreaUp/(area)*ratio;
        (it->second).Maxno = neededno;
        totinc+=neededno;
        
        std::cout<<" We will try to generate  "<<neededno<<"  "<<type<<" protein \n" ;

    }
    //==== end of reading inout file (str) to find number of proteins

    int id=0;
    int s=0;
    //=== checking finding random position for proteins without overlapping. overlapping between proteins not with lipids
    while (totcreated<totinc && s<(pPointUp.size()))
    {
        s++;
        bool accept = true;
        int RNG1=(rand()%totinc)+1;
        int RNG2=(rand()%pPointUp.size());
        int pointid = (pPointUp.at(RNG2))->GetID();
        int tid=0;
        int l=0;
        bool chosen =false;
        int max = 0;
        double R = 0;
        
        
        ///======== here we only try to find one of the protein type randomly based on RNG1
        for ( std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin(); it != m_map_IncID2ProteinLists.end(); it++ )
        {
             max =  (it->second).Maxno;
            int nocreated = (it->second).created;
            //if(l<max && RNG1<=max)
            if(RNG1>l && RNG1<=max+l && nocreated<max)
            {
                tid = (it->second).ID;
                chosen = true;
                
                std::string type = (it->second).ProteinName;
                double area = (m_map_MolName2MoleculesType.at(type)).molarea;
                
                R=sqrt(area/acos(-1));
                
            }
            l=l+max;

        }

        
        
        //===================
        
        for ( std::vector<ExcludedVolumeBeads>::iterator it = m_ExcludeBeads.begin(); it != m_ExcludeBeads.end(); it++ )
        {
            
            Vec3D XP1 = it->X;
            double R1 = it->R;
            
            Vec3D XP2 =(pPointUp.at(RNG2))->GetPos();
            
            if(XP2.dot((XP2-XP1),(XP2-XP1))<(R1+R)*(R1+R))
                accept = false;

            
        }
        
        //======================
        
        if(accept==true)
        {

            double d1= double(rand()%1000)/1000;
            double d2= double(rand()%1000)/1000;
            
            
            //
            Vec3D D(d1,d2,0);
            D=D*(1/(D.norm()));

            Vec3D N = (pPointUp.at(RNG2))->GetNormal();
            Vec3D T1 =   (pPointUp.at(RNG2))->GetP1();
            Vec3D T2 =   (pPointUp.at(RNG2))->GetP2();
            Tensor2 LG = TransferMatLG(N,T1,T2);
            D=LG*D;
            //
            id++;
        inclusion inc(id, tid,pointid,D);
        m_RandomInc.push_back(inc);
            totcreated++;
            (m_map_IncID2ProteinLists.at(tid)).created = (m_map_IncID2ProteinLists.at(tid)).created +1;
            
            
            ExcludedVolumeBeads Ex;
            Ex.X =(pPointUp.at(RNG2))->GetPos();
            Ex.R = R;
            m_ExcludeBeads.push_back(Ex);
            
            
        }
    }
}
bool BackMap::FindProteinList(std::string filename)
{
    //===  this function was updated in Dec, 2023.
    std::cout<<"---> reading protein information from the str file "<<"\n";
    Nfunction f;
    std::ifstream strfile;
    strfile.open(filename.c_str());
    if(strfile.good()==false)
    {
        std::cout<<"---> error (3R22-D): while opening the str file with  name: "<<filename<<"  ; check if the file exist "<<"\n";
        exit(0);
    }
    //== Read the file until reaching the protein section
    bool proteinflag=false;
    std::string str;
    while (true)
    {
        std::getline (strfile,str);
        if(strfile.eof())
        {
            return false;
        }
        str = f.trim(str);
        str.erase(std::remove(str.begin(), str.end(), '['), str.end());
        str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
        str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
        if(str=="ProteinList")
        {
            proteinflag = true;
            break;
        }
    }
    //== now read the protein information
    if(proteinflag==true)
    {
      while (true)
      {
        std::getline (strfile,str);
        std::vector<std::string> Line = f.split(str);
        if(Line.size()!=0)
        {
            if(Line[0]=="End" || Line[0]=="END")
            {
                break;
            }
            else if(Line.size()>=6)
            {
                ProteinList CP;
                CP.ProteinName = Line.at(0);
                CP.ID = f.String_to_Int(Line.at(1));
                CP.Ratio = f.String_to_Double(Line.at(2));
                CP.Phi=f.String_to_Double(Line.at(3));
                CP.Theta=f.String_to_Double(Line.at(4));
                CP.Z0=f.String_to_Double(Line.at(5));
                m_map_IncID2ProteinLists.insert(std::pair<int,ProteinList>(CP.ID, CP));
            }
            else
            {
            std::cout<<"---> error (3R22-E): protein information in the str file is incomplete "<<"\n";
            std::cout<<"---> other possible causs of this error. 1) the protein section in str file does not end with word End "<<"\n";
            exit(0);
            }
        }
      }//while (true)
    }
    else
    {
    std::cout<<"---> warnning: the str file does not contain any information about proteins "<<"\n";
    }
    strfile.close();

    
    if(m_map_IncID2ProteinLists.size()!=0)
    {
      std::cout<<" ---> Protein List and ID have been read from the input file \n";
      for ( std::map<int,ProteinList>::iterator it = m_map_IncID2ProteinLists.begin(); it != m_map_IncID2ProteinLists.end(); it++ )
      {
          std::cout <<"----> "<< it->first  <<" ------> "<< (it->second).ProteinName<< std::endl ;
      }
      std::cout<<"************************************************************** \n";
    }
    return true;
}
Tensor2  BackMap::TransferMatLG(Vec3D Normal, Vec3D t1, Vec3D t2)
{
    Tensor2  GL(t1,t2,Normal);
    Tensor2 LG=GL.Transpose(GL);
    return LG;
    
}

//=== since 2024
// a function that use up the exclsuions by making the area of that specific points zero.
void BackMap::ExcludePointsUsingExclusion(std::vector<exclusion*> &m_pExc, std::vector<point*> &m_pPointUp, std::vector<point*> &m_pPointDown)
{
    
    if(m_pExc.size()!=0)
    {
        std::cout<<" Note: we are excluding points based on exclusion, If it is slow, contact the developer \n";
        
        for ( std::vector<exclusion*>::iterator it = m_pExc.begin(); it != m_pExc.end(); it++ )
        {
            int pointid=(*it)->GetPointID();
            if(pointid<0 || pointid>m_pPointUp.size())
            std::cout<<"error 23456 \n";
            
            point *Up_p1=m_pPointUp.at(pointid);
            Vec3D Pos = Up_p1->GetPos();
            Vec3D N = Up_p1->GetNormal();
            double R = (*it)->GetRadius();
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
                {
                    (*it1)->UpdateArea(0);

                }

                
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
                {
                    (*it1)->UpdateArea(0);
                    
                }
                
            }
        }
    }
}
//=== check if all the incs have been mapped to a proteinlist and the protein list exist in the moltype
bool BackMap::CheckProteinInfo (std::map<int , ProteinList>& plist, std::map<std::string , MolType>& moltype, std::vector<inclusion*> & incs)
{
    
    //== first find all the inc type id and store them
    std::vector<int> inc_typeid;
    for ( std::vector<inclusion*>::iterator it = incs.begin(); it != incs.end(); it++ )
    {
        int id = (*it)->GetTypeID();
        if(count(inc_typeid.begin(), inc_typeid.end(), id)<=0)
            inc_typeid.push_back(id);
    }
    //== check if for each id in incs, there is one in proteinlist
    for ( std::vector<int>::iterator it = inc_typeid.begin(); it != inc_typeid.end(); it++ )
    {
        
        if(plist.find(*it)!=plist.end())
        {
            // now we get the coresponding protein name and check if in the mol def file (either lib or gro files) such a protein exit
            std::string pname = ((plist.find(*it))->second).ProteinName;
            if(moltype.find(pname)!=moltype.end())
            {
                // do we need to check anything more?? for now no
            }
            else
            {
                std::cout<<"error 2522024_a--> there is an inclusion with type id "<<*it<<" and protein name "<<pname<<" however the stracture does not exist, perhaps include a gro file of it \n";
                return false;
            }
        }
        else
        {
            std::cout<<"error 2522024_b--> there is an inclusion with type id "<<*it<<" however it is not mapped to a protein in the str file \n";
            return false;
        }
            
    }
    
    std::cout<<"---> protein information in the str file, lib/gro files and inclusion info matches well \n";
    return true;
}

bool BackMap::PlaceProteins(std::vector<point*> &pPointUp)
{
    for ( std::map<int,ProteinList>::iterator it1 = m_map_IncID2ProteinLists.begin(); it1 != m_map_IncID2ProteinLists.end(); it1++ )
    {
        int plistid=it1->first;
    for ( std::vector<inclusion*>::iterator it = m_pInc.begin(); it != m_pInc.end(); it++ )
    {
        int id = (*it)->GetTypeID();
        if(plistid==id)
        {
            std::string ptype=(m_map_IncID2ProteinLists.at(id)).ProteinName;
            int pointid = (*it)->GetPointID();
            Vec3D  Dir =  (*it)->GetDirection();
            point *Up_p1=pPointUp.at(pointid);
            Vec3D N = Up_p1->GetNormal();
            Vec3D Pos = Up_p1->GetPos();
            Vec3D T1 =   Up_p1->GetP1();
            Vec3D T2 =   Up_p1->GetP2();

            if (m_map_MolName2MoleculesType.count(ptype) == 0)
            {
                std::cout << "---> error: molecule name " <<ptype<<" does not exist in the attached gro files \n";
                return false;
            }
            
            GenProtein(m_map_MolName2MoleculesType.at(ptype), id, Pos, N, Dir, T1,T2);

        }
    }
    plistid++;
    }
    
    return true;
}
bool BackMap::RemovePointsCloseToBeadList(std::vector<point*> &pPointUp, std::vector<point*> &pPointDown, std::vector<bead*> vpbeads, double RCutOff, Vec3D* m_pBox)
{
    
    GenerateUnitCells GCNT(vpbeads, m_pBox,RCutOff,1.0);
    GCNT.Generate();
    // Here, we try to remove the points that are covered by the proteins. We do it by setting the A=0
  
    for ( std::vector<point*>::iterator it = pPointUp.begin(); it != pPointUp.end(); it++ )
    {
        bool rem = false;
        Vec3D Pos1 = (*it)->GetPos();
        Vec3D N =   (*it)->GetNormal();
        Vec3D Pos2 = Pos1 - N*1.5;
        rem = GCNT.anythingaround(Pos1);  //121945
        if(rem==false)
            rem = GCNT.anythingaround(Pos2);
        if(rem==true)
            (*it)->UpdateArea(0);

    }
    if(m_monolayer == false)
    for ( std::vector<point*>::iterator it = pPointDown.begin(); it != pPointDown.end(); it++ )
    {
        bool rem = false;
        Vec3D Pos1 = (*it)->GetPos();
        Vec3D N =   (*it)->GetNormal();
        Vec3D Pos2 = Pos1 - N*1.5;
        rem = GCNT.anythingaround(Pos1);  //121945
        if(rem==false)
            rem = GCNT.anythingaround(Pos2);
        if(rem==true)
            (*it)->UpdateArea(0);

    }
    
    return true;
}


#endif



