#if !defined(AFX_ReadDTSFolder_H_555B2143_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_ReadDTSFolder_H_555B2143_C13C_5648_BF23_124095086234__INCLUDED_


#include "Def.h"
#include "Vec3D.h"
#include "point.h"
#include "inclusion.h"
#include "exclusion.h"


class ReadDTSFolder
{
public:
    
    ReadDTSFolder();
	 ~ReadDTSFolder();


        inline std::vector<point*>  GetUpperPoints()         {return m_pOuterPoint;}
        inline std::vector<point*>  GetInnerPoints()         {return m_pInnerPoint;}
        inline std::vector<inclusion*>  GetInclusion()         {return m_pInclusion;}
        inline std::vector<exclusion*>  GetExclusion()         {return m_pExclusion;}

        inline Vec3D *GetBox()         {return m_pBox;}



public:
    void Read(std::string gmxfilename);




private:


    std::vector<point>  m_OuterPoint;
    std::vector<point>  m_InnerPoint;
    std::vector<inclusion>  m_Inclusion;
    std::vector<exclusion>  m_Exclusion;

    std::vector<point*>  m_pOuterPoint;
    std::vector<point*>  m_pInnerPoint;
    std::vector<inclusion*>  m_pInclusion;
    std::vector<exclusion*>  m_pExclusion;

    Vec3D m_Box;
    Vec3D *m_pBox;


private:
    bool FileExist (const std::string &name);
    std::vector<point> ReadPointObjects(std::string file,int);
    std::vector<inclusion> ReadInclusionObjects(std::string file);
    std::vector<exclusion> ReadExclusionObjects(std::string file);

 //   void ReadPoint(std::string file);
};


#endif
