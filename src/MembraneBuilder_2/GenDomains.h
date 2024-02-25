#if !defined(AFX_GenDomains_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_GenDomains_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <list>
#include <map>
#include <iomanip>
#include <valarray>
#include "Domain.h"
#include "LMatrix.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "point.h"



class GenDomains
{
public:
    
	GenDomains(std::string file, std::vector<point*> point1, std::vector<point*> point2, bool renorm);
	~GenDomains();
    
    inline  std::vector<Domain*> GetDomains()                const  {return m_pAllDomains;}

public:
    void Configure();
private:
    
    
    

private:
    std::vector<point*> m_pPointUp;
    std::vector<point*> m_pPointDown;
    std::vector<Domain*> m_pAllDomains;
    std::vector<Domain> m_AllDomains;

    bool m_Health;



};


#endif
