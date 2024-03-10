#if !defined(AFX_point_H_999B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_point_H_999B21B8_C13C_5648_BF23_124095086234__INCLUDED_


#include "Vec3D.h"
#include "UnitCell.h"

class point
{
public:
    
	point(int id, Vec3D x,std::vector <double> c );
    ~point();


        inline const int GetID()                const  {return m_ID;}
        inline Vec3D  GetPos()                        {return m_Pos;}
        inline std::vector <double>  GetCurvature()                        {return m_C;}
        inline UnitCell *GetpointUnitCell()        const        {return m_PointUnitCell;}




public:
    

  void UpdateBox(Vec3D* z);
  void UpdatePos( Vec3D X);
  void UpdatepointUnitCell(UnitCell * z);

public:


private:

    bool m_UpperLayer;
    Vec3D m_Pos;
    std::vector <double> m_C;
    UnitCell *m_PointUnitCell;
    Vec3D *m_pBox;
    int m_ID;
};


#endif
