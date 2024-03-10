

#include <stdio.h>
#include "GroFile.h"
#include "Nfunction.h"

GroFile::GroFile(std::string gmxfilename)
{

    m_GroFileName = gmxfilename;
    ReadGroFile(m_GroFileName);

}


GroFile::~GroFile()
{
    
}

void GroFile::AddBead(bead b)
{

    m_AllBeads.push_back(b);
}
void GroFile::ReadGroFile(std::string file)
{
    Nfunction f;
    if(file.size()<4)
    {
        file=file+".gro";
    }
    else if(file.at(file.size()-1)=='o' && file.at(file.size()-2)=='r' && file.at(file.size()-3)=='g')
    {
        
    }
    else
    {
        file=file+".gro";
    }
    
    std::string str;
    std::ifstream FGRO;
    FGRO.open(file.c_str());
    std::string title;
    getline (FGRO,title);
    m_Title = title;
    int atomno;
    FGRO>>atomno;
    getline (FGRO,str);
    std::string resname,aname;
    for (int i=0; i<atomno; i++) //NoBeads
    {
        char *cr  = new char [5];
        FGRO.read(cr,5);
        int resid = atoi(cr);
        FGRO.read(cr,5);
        resname = cr;
        resname.erase(remove_if(resname.begin(), resname.end(), isspace), resname.end());
        FGRO.read(cr,5);
        aname = cr;
        aname.erase(remove_if(aname.begin(), aname.end(), isspace), aname.end());
        FGRO.read(cr,5);
        int atomno = atoi(cr);
        char *Xr  = new char [8];
        FGRO.read(Xr,8);
        double X=atof(Xr);
        FGRO.read(Xr,8);
        double Y=atof(Xr);
        FGRO.read(Xr,8);
        double Z=atof(Xr);
        getline (FGRO,title);
        
        bead Be(i, aname, aname, resname, resid, X, Y, Z);
        m_AllBeads.push_back(Be);
        
    }
    float Lx,Ly,Lz;
    char *Lr  = new char [10];
    FGRO.read(Lr,10);
    Lx=atof(Lr);
    FGRO.read(Lr,10);
    Ly=atof(Lr);
    FGRO.read(Lr,10);
    Lz=atof(Lr);
    FGRO.close();
    
    m_Box(0)=Lx; m_Box(1)=Ly; m_Box(2)=Lz;
    m_pBox = &m_Box;
    double xcm =0;
    double ycm =0;
    double zcm =0;

    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        (*it).UpdateBox(m_pBox);
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        
        xcm+=x/double(m_AllBeads.size());
        ycm+=y/double(m_AllBeads.size());
        zcm+=z/double(m_AllBeads.size());


        
    }
    
   /* for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        (*it).UpdateBox(m_pBox);
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos(); 
        (*it).UpdateXPos(x-xcm);
        (*it).UpdateYPos(y-ycm);
        (*it).UpdateZPos(z-zcm);
        
        
    }*/
    
    
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        m_pAllBeads.push_back(&(*it));
  
    }
    
}

void GroFile::WriteGroFile(std::string file)
{

    if(file.size()<4)
    {
        file=file+".gro";
    }
    else if(file.at(file.size()-1)=='o' && file.at(file.size()-2)=='r' && file.at(file.size()-3)=='g')
    {
        
    }
    else
    {
        file=file+".gro";
    }

    
    FILE *fgro;
    fgro = fopen(file.c_str(), "w");
    
    
    /// resid  res name   noatom   x   y   z
    const char* Title="dmc gmx file handler";
    int Size=m_AllBeads.size();
    
    fprintf(fgro,  "%s\n",Title);
    fprintf(fgro, "%5d\n",Size);
    int i=0;
    for (std::vector<bead>::iterator it = m_AllBeads.begin() ; it != m_AllBeads.end(); ++it)
    {
        
        i++;
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        
        
        const char* A1=((*it).GetResName()).c_str();
        const char* A2=((*it).GetBeadName()).c_str();
        int resid=(*it).GetResid();
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid,A1,A2,i,x,y,z );

    }
    

    fprintf(fgro,  "%10.5f%10.5f%10.5f\n",m_Box(0),m_Box(1),m_Box(2) );
    fclose(fgro);
}


