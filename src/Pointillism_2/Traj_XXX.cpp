
#include <time.h>
#include <iomanip>
#include "Traj_XXX.h"


Traj_XXX::Traj_XXX(Vec3D *pBox)
{
 m_pBox=pBox;
m_Condition=true;
    Nfunction f;
    m_tsiPrecision = G_tsiPrecision ; //f.Int_to_String(18)+"."+f.Int_to_String(10);
}

Traj_XXX::~Traj_XXX()
{
    
}
void Traj_XXX::WriteTSI(int step ,  std::string filename , std::vector< vertex* > pver, std::vector< triangle* > ptriangle,  std::vector< inclusion* > pinc , std::vector< exclusion* > pexc)
{
    FILE * output;
    output = fopen((filename).c_str(), "w");
    std::string format = "%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    const char* version=SoftWareVersion;
    fprintf(output,"%s\n",version);
//------
    const char* box="box";
    format = "%s%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    fprintf(output,format.c_str(),box,(*m_pBox)(0),(*m_pBox)(1),(*m_pBox)(2));

    const char* ver="vertex";
    int size=pver.size();
    fprintf(output,"%s%20d\n",ver,size);
    format = "%5d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<vertex *>::iterator it = pver.begin() ; it != pver.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetVID(),(*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
    
    const char* tri="triangle";
    size = ptriangle.size();
    fprintf(output,"%s%20d\n",tri,size);
    for (std::vector<triangle *>::iterator it = ptriangle.begin() ; it != ptriangle.end(); ++it)
        fprintf(output,"%10d%10d%10d%10d\n",(*it)->GetTriID(),((*it)->GetV1())->GetVID(),((*it)->GetV2())->GetVID(),((*it)->GetV3())->GetVID());
    
    
    const char* inc="inclusion";
    size = pinc.size();
    fprintf(output,"%s%20d\n",inc,size);
    format = "%10d%10d%10d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    

    
    for (std::vector<inclusion *>::iterator it = pinc.begin() ; it != pinc.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetID(),((*it)->GetInclusionTypeID()),((*it)->Getvertex())->GetVID(),((*it)->GetLDirection())(0),((*it)->GetLDirection())(1));
    
    if(pexc.size()!=0)
    {
        const char* exc="exclusion";
        size = pexc.size();
        fprintf(output,"%s%20d\n",exc,size);
        format = "%10d%10d%"+m_tsiPrecision+"lf\n";
        for (std::vector<exclusion *>::iterator it = pexc.begin() ; it != pexc.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetID(),((*it)->Getvertex())->GetVID(),(*it)->GetRadius());
    }
    
    
    
    
    fclose(output);
    
}

