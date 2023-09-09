

#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdlib>
#include "Edit_configuration.h"
#include "WriteFiles.h"
#include "help.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "VMDOutput.h"
#include "Trajectory.h"
#include "Topology.h"
#include "Surface_Mosaicing.h"
#include "VertexMove.h"
#include "Traj_XXX.h"
#include "In_OR_Out.h"


/*
 This has been updated in Sept 2023. This software is now both backmapping and also input generator.
 
 this class is used to edit output configuration from a DTS simulations. Mosaicing increases this number of the verteices to be used for backmapping,
then generate two surfaces that can be used to create a bilayer.
 
 first rescale the trjectory file,
 Mosaic the surface few times.
 Then generate two surfaces .
 
 */
Edit_configuration::Edit_configuration( std::vector <std::string> Arguments)
{

    // Initialize the Variables to their default values
    InitializeVariables();
    // updates variables based on the command line arguments
    UpdateVariables(Arguments);

    // generating the mesh
    std::string domyfile;
    CreateMashBluePrint BluePrint;
    MeshBluePrint meshblueprint;
    meshblueprint = BluePrint.MashBluePrintFromInput_Top(domyfile,m_MeshFileName);
    std::vector<double> x;
    //m_Mesh.GenerateMesh(meshblueprint,0,0,x);
    m_pMesh = &m_Mesh;
    m_pCurvatureCalculations = &m_CurvatureCalculations;
    // Mesh is ready
    

        m_pBox=m_pMesh->m_pBox;
    

//// do the jobs
  if(m_TaskName == "in_out")
  {
      In_OR_Out in_out(m_MeshFileName);
  }
  else if(m_TaskName=="check")
  {
        check(m_MeshFileName);
  }
  else if(m_TaskName=="minimize")
  {
        Minimize(m_MeshFileName);
  }
  else if(m_TaskName=="vertexinfo")
  {
          // a new feature in TS2CG1.2 and above. Finding geometric info on a vertex
          VertexInfo(m_MeshFileName);
  }
  else if (m_TaskName=="PLM")
  {
        double H=m_BilayerThickness/2.0;
        const int dir_err = system(("mkdir -p "+m_Folder).c_str());
        if (-1 == dir_err)
        {
            std::cout<<"error--> creating directory  "<<m_Folder<<"\n";
            exit(1);
        }
        const int dir_err2 = system(("mkdir -p "+m_Folder+"visualization_data").c_str());
        if (-1 == dir_err2)
        {
            std::cout<<"error--> creating directory  visualization_data"<<"\n";
            exit(1);
        }
            BackMapOneLayer(1 , m_MeshFileName, H);
            if(m_monolayer==0)
            BackMapOneLayer(-1 , m_MeshFileName, H);
  }
  else
    std::cout<<" error--> unrecognized Task \n";
}
Edit_configuration::~Edit_configuration()
{
    
}

void Edit_configuration::Rescaling(Vec3D zoom )
{
    Vec3D L((*m_pBox)(0),(*m_pBox)(1),(*m_pBox)(2));
    (*m_pBox)(0) = zoom(0)*L(0);
    (*m_pBox)(1) = zoom(1)*L(1);
    (*m_pBox)(2) = zoom(2)*L(2);
    
    for (std::vector<vertex *>::iterator it = (m_pMesh->m_pActiveV).begin() ; it != (m_pMesh->m_pActiveV).end(); ++it)
    {
        double x = (*it)->GetVXPos();
        double y = (*it)->GetVYPos();
        double z = (*it)->GetVZPos();
        (*it)->UpdateVXPos(zoom(0)*x);
        (*it)->UpdateVYPos(zoom(1)*y);
        (*it)->UpdateVZPos(zoom(2)*z);
        
    }
    return;
}
bool Edit_configuration::FileExist (const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}
void  Edit_configuration::UpdateGeometry()
{

    
    for (std::vector<triangle *>::iterator it = (m_pMesh->m_pActiveT).begin() ; it != (m_pMesh->m_pActiveT).end(); ++it)
    (*it)->UpdateNormal_Area(m_pBox);

    //===== Prepare links:  normal vector and shape operator
    for (std::vector<links *>::iterator it = (m_pMesh->m_pHL).begin() ; it != (m_pMesh->m_pHL).end(); ++it)
    {
            (*it)->UpdateNormal();
            (*it)->UpdateShapeOperator(m_pBox);
    }

    //======= Prepare vertex:  area and normal vector and curvature of surface vertices not the edge one
    for (std::vector<vertex *>::iterator it = (m_pMesh->m_pSurfV).begin() ; it != (m_pMesh->m_pSurfV).end(); ++it)
        m_pCurvatureCalculations->SurfVertexCurvature(*it);
        
    //====== edge links should be updated
    for (std::vector<links *>::iterator it = (m_pMesh->m_pEdgeL).begin() ; it != (m_pMesh->m_pEdgeL).end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = (m_pMesh->m_pEdgeV).begin() ; it != (m_pMesh->m_pEdgeV).end(); ++it)
            m_pCurvatureCalculations->EdgeVertexCurvature(*it);
    

    
}
bool Edit_configuration::check(std::string file){

    UpdateGeometry();

        double Tarea  = 0;
        double totalgaussianC = 0;
        double pi = acos(-1);
        for (std::vector<vertex *>::iterator it = (m_pMesh->m_pSurfV).begin() ; it != (m_pMesh->m_pSurfV).end(); ++it)
        {
            
            Vec3D X((*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
            Vec3D normal = (*it)->GetNormalVector();
            Tarea+= (*it)->GetArea();
            std::vector <double> C = (*it)->GetCurvature();
            totalgaussianC+= (C.at(0))*(C.at(1))*((*it)->GetArea());
        }
        int nv = (m_pMesh->m_pSurfV).size();
        int nt = (m_pMesh->m_pActiveT).size();
        int nl = (m_pMesh->m_pHL).size();
        std::cout<<" ************************************************************************************** \n";
        std::cout<<" Total area of the surface is "<<Tarea<<" nm^2 you can use -rescalefactor to increase it \n";
        //std::cout<<" total gaussian curvature /2PI is: "<<totalgaussianC/(2*pi)<<" \n";
       // std::cout<<" total mean curvature /2PI is: "<<totalmeanC/(2*pi)<<" \n";
        std::cout<<" Euler characteristic: "<<nv+nt-nl<<" \n";
        std::cout<<" topology genus "<<(2-(nv+nt-nl))/2<<" \n";


    return true;
}
void Edit_configuration::VertexInfo(std::string file){

        UpdateGeometry();
        int index;
        std::cout<<" Please provide the index of the vertex (note, it will start from zero): \n";
        std::cin>>index;
        
        if(index<0 || index>(m_pMesh->m_pActiveV).size())
        {
            std::cout<<" ---> error:the index number is invalid \n";
            exit(0);
        }
        vertex *pv = (m_pMesh->m_pActiveV).at(index);
        

        double x = pv->GetVXPos();
        double y = pv->GetVYPos();
        double z = pv->GetVZPos();
        double A = pv->GetArea();
        Tensor2  L2GT = pv->GetL2GTransferMatrix();
        Tensor2  G2LT = pv->GetG2LTransferMatrix();
        Vec3D normal  = pv->GetNormalVector();
        std::vector <double> C = pv->GetCurvature();
        
        std::cout<<" vertex position is "<<x<<"  "<<y<<"  "<<z<<" \n ";
        std::cout<<" vertex area is "<<A<<" \n ";
        std::cout<<" vertex normal is "<<normal(0)<<"  "<<normal(1)<<"  "<<normal(2)<<" \n ";

        Vec3D t1(1,0,0);
        Vec3D t2(0,1,0);
        Vec3D n(0,0,1);
        
        Vec3D T1 = L2GT*t1;
        Vec3D T2 = L2GT*t2;
        Vec3D N = L2GT*n;
        std::cout<<" first direction "<<T1(0)<<"  "<<T1(1)<<"  "<<T1(2)<<" \n ";
        std::cout<<" second direction "<<T2(0)<<"  "<<T2(1)<<"  "<<T2(2)<<" \n ";
        std::cout<<" normal "<<N(0)<<"  "<<N(1)<<"  "<<N(2)<<" \n ";
}


/// ====== Updates in 2023
void Edit_configuration::UpdateVariables(std::vector <std::string> &Arguments)
{
    Nfunction f;    // an object to use some of pre-made functions
    std::ofstream log;
    log.open("plm.log");
    log<<" outputs have been generated by below command \n";

    for (long i=0;i<Arguments.size();i++)
    {
        log<<Arguments[i]<<"  ";
    }
// read the arguments in the command line and update the variables
    for (int i=1;i<Arguments.size();i=i+2)
    {
        if(Arguments[i]==Def_TSfile)
        {
            m_MeshFileName=Arguments[i+1];    // ts file name, *.q, *.tsi, *.dat
            
        }
        else if(Arguments[i]==Def_HelpCall)
        {
            help helpmessage(SoftWareVersion, Arguments[0]);
            exit(0);
        }
        else if(Arguments[i]==Def_BilayerThickness)
        {
            m_BilayerThickness=f.String_to_Double(Arguments.at(i+1));  //bilayer thickness
        }
        else if(Arguments[i]==Def_AlgType)
        {
            m_MosAlType=(Arguments.at(i+1));  // algorithm type, Type1 and Type2
        }
        else if(Arguments[i]==Def_rescalefactor)
        {
            m_Zoom(0)=f.String_to_Double(Arguments.at(i+1));
            m_Zoom(1)=f.String_to_Double(Arguments.at(i+2));
            m_Zoom(2)=f.String_to_Double(Arguments.at(i+3));
            i=i+2;
        }
        else if(Arguments[i]==Def_AreaPerLipid)
        {
            m_AP=f.String_to_Double(Arguments[i+1]);
        }
        else if(Arguments[i]==Def_DegreeOfMeshing)
        {
            m_Iteration=f.String_to_Int(Arguments.at(i+1));
        }
        else if(Arguments[i]==Def_TaskName) //"-r"
        {
            m_TaskName=Arguments[i+1];
        }
        else if(Arguments[i]==Def_OutputFolderName) // "-o"
        {
            m_Folder=Arguments[i+1];
        }
        else if(Arguments[i]==Def_SmoothingFlag)// "-smooth"
        {
            m_smooth=true;
            i=i-1;
        }
        else if(Arguments[i]==Def_resizebox)// "-resizebox"
        {
            m_FindnewBox=true;
            i=i-1;
        }
        else if(Arguments[i]==Def_Monolayer)   // "-monolayer"
        {
            m_monolayer=f.String_to_Int(Arguments[i+1]);
            m_BilayerThickness = 0;
        }
        else
        {
            std::cout<<"error---> unrecognized argument || "<<Arguments[i]<<" || \n";
            log<<"error---> unrecognized argument || "<<Arguments[i]<<" || \n";
            exit(0);
        }
    }
    
    if(f.FileExist(m_MeshFileName)==false)
    {
    std::cout<<"error--> ts file with the name "<<m_MeshFileName<<" does not exist in the folder \n";
    log<<"error--> ts file with the name "<<m_MeshFileName<<" does not exist in the folder \n";
        exit(0);
    }
    
    log.close();
    return;
}
void Edit_configuration::InitializeVariables()
{
    //== Initialize  variables
        m_smooth = false;
        m_Folder = "point";
        m_TaskName = "PLM";
        m_monolayer = 0;
        m_FindnewBox = false; // if true, find a box as small as possible
        m_MosAlType = "Type1"; // algorithm type, Type1 and Type2
        m_BilayerThickness = 3.8;
        m_Zoom(0) = 1;
        m_Zoom(1) = 1;
        m_Zoom(2) = 1;
        m_Iteration = -1;
        m_MeshFileName = "TS.q";
        m_AP = 0.62;
    
    
    return;
}
double  Edit_configuration::PPBCM_Cluster(double Lx, std::vector <double> X)
{
    double cm=0;
    
    
    std::vector <double> X1,X2;
    for (std::vector<double>::iterator it = X.begin() ; it != X.end(); ++it)
    {
        if((*it)>Lx/2)
        {
            X2.push_back((*it));
        }
        else
        {
            X1.push_back((*it));
        }
        
    }
    double cm1=0;
    double cm2=0;
    for (std::vector<double>::iterator it = X1.begin() ; it != X1.end(); ++it)
    {
        cm1+=(*it)/X1.size();
        
    }
    for (std::vector<double>::iterator it = X2.begin() ; it != X2.end(); ++it)
    {
        cm2+=(*it)/X2.size();
        
    }
    
        bool boxout=false;
    
    if(fabs(cm2-cm1)>Lx/2)
        boxout=true;
    if(boxout==false)
    {
        cm=(cm2*X2.size()+cm1*X1.size())/X.size();
    }
    else
    {
        cm=(cm2*X2.size()+(cm1+Lx)*X1.size())/X.size();
    }

    return cm;
    
}
void Edit_configuration::UpdateBoxSize()
{
    // Finding the center of geometry
    Vec3D CM(0,0,0);
    int totalvNo = (m_pMesh->m_pActiveV).size();
    for (std::vector<vertex *>::iterator it1 = (m_pMesh->m_pActiveV).begin() ; it1 != (m_pMesh->m_pActiveV).end(); ++it1)
    {
        Vec3D X ((*it1)->GetVXPos(),(*it1)->GetVYPos(),(*it1)->GetVZPos());
        CM = CM+X*(1.0/double(totalvNo));
    }
    // Finding the max distance from the CMG
    Vec3D MaxB(0,0,0);
    for (std::vector<vertex *>::iterator it = (m_pMesh->m_pActiveV).begin() ; it != (m_pMesh->m_pActiveV).end(); ++it)
    {

            double dx = (*it)->GetVXPos()-CM(0);
            double dy = (*it)->GetVYPos()-CM(1);
            double dz = (*it)->GetVZPos()-CM(2);
        
            if(MaxB(0)<fabs(dx))
                MaxB(0)=fabs(dx);
            if(MaxB(1)<fabs(dy))
                MaxB(1)=fabs(dy);
            if(MaxB(2)<fabs(dz))
                MaxB(2)=fabs(dz);

    }
    double H = m_BilayerThickness/2;
    Vec3D DB ((6+H)/m_Zoom(0),(6+H)/m_Zoom(1),(6+H)/m_Zoom(2));
    (*m_pBox) = MaxB*2+DB;
    for (std::vector<vertex *>::iterator it1 = (m_pMesh->m_pActiveV).begin() ; it1 != (m_pMesh->m_pActiveV).end(); ++it1)
    {
        Vec3D X ((*it1)->GetVXPos(),(*it1)->GetVYPos(),(*it1)->GetVZPos());
        (*it1)->UpdateVXPos((*it1)->GetVXPos()-CM(0)+((*m_pBox)(0))*0.5);
        (*it1)->UpdateVYPos((*it1)->GetVYPos()-CM(1)+((*m_pBox)(1))*0.5);
        (*it1)->UpdateVZPos((*it1)->GetVZPos()-CM(2)+((*m_pBox)(2))*0.5);

    }
    
}
//=== the backmapping function
//==================================
void Edit_configuration::BackMapOneLayer(int layer , std::string file, double H)
{
   
    

  if(m_FindnewBox==true)
      UpdateBoxSize();
    
   UpdateGeometry();
    
    // this means that we need to find the optiomal iteration value
    double totalvNo = double((m_pMesh->m_pActiveV).size());
    if(m_Iteration==-1 )
    {
            double LargeZoom = m_Zoom(0);
            if(LargeZoom<m_Zoom(1))
            LargeZoom = m_Zoom(1);
            if(LargeZoom<m_Zoom(2))
            LargeZoom = m_Zoom(2);
        
            double Tarea  = 0;
            for (std::vector<vertex *>::iterator it = (m_pMesh->m_pActiveV).begin() ; it != (m_pMesh->m_pActiveV).end(); ++it)
            Tarea+= (*it)->GetArea();
            double requiredNo = Tarea*LargeZoom*LargeZoom/(m_AP*totalvNo);
            m_Iteration = int (log(requiredNo)/log(4))+1;
            if(m_Iteration<5)
            std::cout<<"We will increase the number of the available points by 4^"<<m_Iteration<<"\n";
            else
            std::cout<<"The requested rescaling requires  4^"<<m_Iteration<<" points, this may take some time to finish \n";
    }

    
    Rescaling(m_Zoom);

    // Moving each vertex in the direction of the normal vector
    for (std::vector<vertex *>::iterator it = (m_pMesh->m_pActiveV).begin() ; it != (m_pMesh->m_pActiveV).end(); ++it)
    {
        double x = (*it)->GetVXPos();
        double y = (*it)->GetVYPos();
        double z = (*it)->GetVZPos();
        Vec3D normal = (*it)->GetNormalVector();
        x=x+layer*H*(normal(0));
        y=y+layer*H*(normal(1));
        z=z+layer*H*(normal(2));
        
        (*it)->UpdateVXPos(x);
        (*it)->UpdateVYPos(y);
        (*it)->UpdateVZPos(z);
    }
    UpdateGeometry();
    std::vector <Surface_Mosaicing> Vmos;
    /* for (int j=0;j<m_Iteration;j++)
    {
        Surface_Mosaicing  MOS(m_MosAlType,m_smooth);
        Vmos.push_back(MOS);

    }
  


    
    






#if BACKMAP == Enabled

    //

    for (int j=0;j<m_Iteration;j++)
    {
        //std::cout<<" we are here 10 \n";
        // Here will cause error when
        (Vmos.at(j)).PerformMosaicing(m_pBox,  m_pAllV ,  m_pAllT , m_pAllLinks, m_pInc,m_Iteration-j);
        //std::cout<<" we are here 11 \n";

        m_pAllV.clear();
        m_pAllT.clear();
        m_pAllLinks.clear();
        m_pHalfLinks1.clear();
        
        m_pAllLinks = (Vmos.at(j)).GetLinks();
        m_pAllV = (Vmos.at(j)).GetVertexs();
        m_pAllT = (Vmos.at(j)).GetTriangles();
        m_pHalfLinks1 = (Vmos.at(j)).GetHLinks();
        
    }
#elif
    std::cout<<"turn on Backmap Flag in the Def "<<std::endl;
#endif
  
    double Lx=(*m_pBox)(0);
    double Ly=(*m_pBox)(1);
    double Lz=(*m_pBox)(2);
    
  
    {
        Vec3D BoxSides=(*m_pBox);
        std::string filename;
        
        if(layer==1)
        filename=m_Folder+"visualization_data/Upper";
        if(layer==-1)
        filename=m_Folder+"visualization_data/Lower";
        VMDOutput GRO(BoxSides, m_pAllV , m_pHalfLinks1, filename);
        GRO.WriteGro();
        GRO.WriteGro2();
        WriteFiles vtu(m_pBox);
        std::string fi;
        if(layer==1)
        fi=m_Folder+"visualization_data/Upper.vtu";
        if(layer==-1)
        fi=m_Folder+"visualization_data/Lower.vtu";
        vtu.Writevtu(m_pAllV,m_pAllT,m_pAllLinks,fi);
        TSI.WriteTSI(0,"extended.tsi",m_pAllV,m_pAllT,m_pInc,m_pExc);
    }
    
    //=============
    
    std::string     UFUpper = m_Folder+"/OuterBM.dat";
    std::string     UFInner = m_Folder+"/InnerBM.dat";
    
    
    
    FILE *BMFile1;
    if(layer==1)
    BMFile1 = fopen(UFUpper.c_str(), "w");
    FILE *BMFile2;
    if(layer==-1)
    BMFile2 = fopen(UFInner.c_str(), "w");
    
    const char* Cbox="Box";
    
    if(layer==1)
    fprintf(BMFile1,  "%s%12.3f%12.3f%12.3f\n",Cbox,Lx,Ly,Lz);
    if(layer==2)
    fprintf(BMFile2,  "%s%12.3f%12.3f%12.3f\n",Cbox,Lx,Ly,Lz);
    
    const char* STR1="< Point NoPoints";
    const char* STR2=">";
    
    int NoPoints=m_pAllV.size();
    if(layer==1)
    fprintf(BMFile1,  "%s%10d%s\n",STR1,NoPoints,STR2);
    if(layer==-1)
    fprintf(BMFile2,  "%s%10d%s\n",STR1,NoPoints,STR2);
    
    const char* Cont="< id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2  >";
    if(layer==1)
    fprintf(BMFile1,  "%s\n",Cont);
    if(layer==-1)
    fprintf(BMFile2,  "%s\n",Cont);
    
    {
        const char* lay="< Outer >";
        if(layer==1)
        fprintf(BMFile1,  "%s\n",lay);
    }
    {
        const char* lay="< Inner >";
        if(layer==-1)
        fprintf(BMFile2,  "%s\n",lay);
    }
    
    
    
    int i=0;
    double dr=1;
    if (m_monolayer==-1)
    dr=-1;
    for (std::vector<vertex *>::iterator it = m_pAllV.begin() ; it != m_pAllV.end(); ++it)
    {
        
        double area=(*it)->GetArea();
        Vec3D normal=(*it)->GetNormalVector();
        normal = normal*(layer)*(dr);
        std::vector <double> curvature = (*it)->GetCurvature();
        //    inclusion* inc = (*it)->GetInclusion();
        Tensor2  L2G = (*it)->GetL2GTransferMatrix();
        Vec3D LD1(1,0,0);
        Vec3D GD1 = L2G*LD1;
        Vec3D LD2(0,1,0);
        Vec3D GD2 = L2G*LD2;
        double x=(*it)->GetVXPos();
        double y=(*it)->GetVYPos();
        double z=(*it)->GetVZPos();
        int domain = (*it)->GetDomainID();

        
        
        if(layer==1) fprintf(BMFile1,"%5d%5d%10.3f%10.3f%10.3f%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",i,domain,area,x,y,z,normal(0),normal(1),normal(2),GD1(0),GD1(1),GD1(2),GD2(0),GD2(1),GD2(2),curvature.at(0),curvature.at(1));
        
        if(layer==-1) fprintf(BMFile2,"%5d%5d%10.3f%10.3f%10.3f%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",i,domain,area,x,y,z,normal(0),normal(1),normal(2),GD1(0),GD1(1),GD1(2),GD2(0),GD2(1),GD2(2),curvature.at(0),curvature.at(1));
        
        i++;
        
    }
    
if(layer==1)
{
    
    
   if (m_pInc.size()!=0){
    FILE *IncFile;
    IncFile = fopen((m_Folder+"/IncData.dat").c_str(), "w");
    
    
    
    const char* CHAR1 ="< Inclusion NoInc   ";
    NoPoints=m_pInc.size();
    const char* CHAR2 ="   >";
    
    fprintf(IncFile,  "%s%5d%s\n",CHAR1,NoPoints,CHAR2);
    
    const char* CHAR3 ="< id typeid pointid lx ly lz  >";
    fprintf(IncFile,  "%s\n",CHAR3);
    
    
    
    i=0;
    
    for (std::vector<inclusion *>::iterator it = m_pInc.begin() ; it != m_pInc.end(); ++it)
    {
        
        int intypeid = ((*it)->GetInclusionType())->ITid;
        vertex* ver = (*it)->Getvertex();
        Tensor2  L2G = ver->GetL2GTransferMatrix();
        Vec3D LD = (*it)->GetLDirection();
        //std::cout<<LD(0)<<" LD "<<LD(1)<<"  "<<LD(2)<<"  \n";
        Vec3D GD = L2G*LD;
        //Tensor2  G2L = L2G.Transpose(L2G);
        //std::cout<<GD(0)<<" GD "<<GD(1)<<"  "<<GD(2)<<"  \n";
        //std::cout<<(G2L*GD)(0)<<" LGD "<<(G2L*GD)(1)<<"  "<<(G2L*GD)(2)<<"  \n";

        int verid=ver->GetVID();
        fprintf(IncFile,  "%12d%12d%12d%8.3f%8.3f%8.3f\n",i,intypeid,verid,GD(0),GD(1),GD(2));
        i++;
    }
   }
    //==== We write Exclusion data
    
    if(m_pExc.size()!=0)
    {
    FILE *ExcFile;
    ExcFile = fopen((m_Folder+"/ExcData.dat").c_str(), "w");
    
    
    
    const char* CHAR1 ="< Exclusion NoExc   ";
    NoPoints=m_pExc.size();
    const char* CHAR2 ="   >";
    
    fprintf(ExcFile,  "%s%5d%s\n",CHAR1,NoPoints,CHAR2);
    
    const char* CHAR3 ="< id  pointid r >";
    fprintf(ExcFile,  "%s\n",CHAR3);
    
    
    
    i=0;
    
    for (std::vector<exclusion *>::iterator it = m_pExc.begin() ; it != m_pExc.end(); ++it)
    {
        
        vertex* ver = (*it)->Getvertex();
        int verid=ver->GetVID();
        double R=(*it)->GetRadius();
        fprintf(ExcFile,  "%5d%10d%8.3f\n",i,verid,R);
        i++;
    }
    }
}
    
    
    */
    
}
void Edit_configuration::Minimize(std::string file){
    
   
//=================================
    
    std::cout<<" We are trying to smoothing the TS for low scale geo \n";
    m_minRoughness = 0.4;

        UpdateGeometry();
        srand (8753);
 /*
        Surface_Mosaicing S(m_MosAlType, true);
        double DR=0.02;
        double DX = 0.1;
        for (std::vector<links *>::iterator it = m_pHalfLinks1.begin() ; it != m_pHalfLinks1.end(); ++it)
        {
            double linklength,midpointdistance;
            S.RoughnessOfALink((*it), &linklength, &midpointdistance);
            
            double oldroughness = midpointdistance-m_minRoughness*linklength;
            if(oldroughness>0)
            {
                vertex * pv1=(*it)->GetV1();
                vertex * pv2=(*it)->GetV2();
                Vec3D P10(pv1->GetVXPos(),pv1->GetVYPos(),pv1->GetVZPos());
                Vec3D P20(pv2->GetVXPos(),pv2->GetVYPos(),pv2->GetVZPos());

                pv1->UpdateOwnInclusion(true);
                pv2->UpdateOwnInclusion(true);

                bool stop=false;
                std::cout<<" Link with ID "<<(*it)->GetID()<<" and ver1 "<<pv1->GetVID()<<"  ver2 "<<pv1->GetVID()<<"\n";
                int loop= 0;
                while(stop==false)
                {
                    loop++;
                    double dx=1-2*((double(rand()%2000000)/2000000.0));           // Inside a cube with the side length of R
                    double dy=1-2*((double(rand()%2000000)/2000000.0));
                    double dz=1-2*((double(rand()%2000000)/2000000.0));
                    double chosever=((double(rand()%2000000)/2000000.0));

                    std::vector <double> C1=pv2->GetCurvature();
                    std::vector <double> C2=pv1->GetCurvature();
                    double A1=pv2->GetArea();
                    double A2=pv1->GetArea();
                    double E1=A1*(C1.at(0)+C1.at(1))*(C1.at(0)+C1.at(1))+A2*(C2.at(0)+C2.at(1))*(C2.at(0)+C2.at(1));


                    vertex *TV;
                    if(chosever>0.5)
                        TV=pv2;
                    else
                        TV=pv1;
                    
                    

                        VertexMove  vmove(TV, DR*dx, DR*dy, DR*dz,m_pBox);
                        vmove.Move();
                        C1.clear();C2.clear();
                        C1=pv2->GetCurvature();
                        C2=pv1->GetCurvature();
                        A1=pv2->GetArea();
                        A2=pv1->GetArea();
                        double E2=A1*(C1.at(0)+C1.at(1))*(C1.at(0)+C1.at(1))+A2*(C2.at(0)+C2.at(1))*(C2.at(0)+C2.at(1));
                        
                        S.RoughnessOfALink((*it), &linklength, &midpointdistance);
                        double roughness = midpointdistance-m_minRoughness*linklength;
                    
                    Vec3D P1(pv1->GetVXPos(),pv1->GetVYPos(),pv1->GetVZPos());
                    Vec3D P2(pv2->GetVXPos(),pv2->GetVYPos(),pv2->GetVZPos());
                    
                    (P1-P10).norm();
                    (P2-P20).norm();
                    
                        if(roughness<oldroughness  && ((P1-P10).norm())<DX && ((P2-P20).norm())<DX && E2-E1<0.1)
                            oldroughness = roughness;
                        else
                            vmove.RejectMove();

                        if((roughness<0 ||loop>10000))
                        {
                            std::cout<<" the roughness is reduced to "<<roughness<<"  with "<<loop<<" iterration \n";
                            stop=true;
                        }




                    
                }

            }

        }
        

        

        
    }
    
    
    WriteFiles vtu(m_pBox);
    std::string fi = "Smooth.vtu";
    vtu.Writevtu(m_pAllV,m_pAllT,m_pAllLinks,fi);
    
    std::ofstream output;
    output.open("smoothtopo.q");
    
    
    output<<std::fixed;
    output<<std::setprecision( Precision )<<(*m_pBox)(0)<<"   "<<(*m_pBox)(1)<<"   "<<(*m_pBox)(2)<<"   \n";
    output<<m_pAllV.size()<<"\n";
    for (int i=0;i<m_pAllV.size();i++)
    {
        vertex* a=m_pAllV.at(i);
        output<<std::setprecision( Precision )<<i<<"  "<<a->GetVXPos()<<" "<<a->GetVYPos()<<" "<<a->GetVZPos()<<" "<<0<<"\n";

    }
    output<< m_pAllT.size()<<"\n";
    for (int i=0;i<m_pAllT.size();i++)
    {

        triangle* a=m_pAllT.at(i);
        if(a->GetRepresentation()==true)
            output<<i<<"   "<<(a->GetV1())->GetVID()<<" "<<(a->GetV2())->GetVID()<<" "<<(a->GetV3())->GetVID()<<" "<<1<<"\n";
        else
            output<<i<<"   "<<(a->GetV1())->GetVID()<<" "<<(a->GetV2())->GetVID()<<" "<<(a->GetV3())->GetVID()<<" "<<0<<"\n";
    }

*/
}
