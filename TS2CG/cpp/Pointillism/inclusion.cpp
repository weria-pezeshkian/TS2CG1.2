 #if !defined(AFX_inclusion_CPP_7F4A21C7_C13Q_8823_BF2E_124095086234__INCLUDED_)
#define AFX_inclusion_CPP_7F4A21C7_C13Q_8823_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include "inclusion.h"
#include "vertex.h"
#include "Nfunction.h"

inclusion::inclusion(int id)
{
m_ID=id;
m_TypeID = 0;
}

inclusion::~inclusion()
{
    
}
void inclusion::UpdateInclusionTypeID(int Typeid)
{
    m_TypeID = Typeid;
}
void inclusion::Updatevertex(vertex * v)
{
    m_pvertex = v;
}
void inclusion::UpdateLocalDirection(Vec3D v)
{
    m_LDirection = v;
}
#endif



