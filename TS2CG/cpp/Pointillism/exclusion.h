#if !defined(AFX_exclusion_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_exclusion_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_
#include "SimDef.h"
#include "Vec3D.h"
#include "vertex.h"
#include "Nfunction.h"

class vertex;
class exclusion
{
public:
    
	exclusion(int id);
	 ~exclusion();


	    inline const int GetID()                const  {return m_ID;}
        inline vertex* Getvertex()              {return m_pvertex;}
        inline double GetRadius()               {return m_R;}


public:
    

  void Updatevertex(vertex * );
  void UpdateRadius(double );

public:


private:



    int m_ID;
    double m_R;
    vertex *m_pvertex;
    









};


#endif
