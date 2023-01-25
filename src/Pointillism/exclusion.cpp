

#include "exclusion.h"


exclusion::exclusion(int id)
{


m_ID=id;
m_R = 0.0;



}

exclusion::~exclusion()
{
    
}

void exclusion::Updatevertex(vertex * v)
{
    m_pvertex = v;
}
void exclusion::UpdateRadius(double r)
{
    m_R = r;
}



