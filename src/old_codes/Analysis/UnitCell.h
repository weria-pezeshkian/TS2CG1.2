#if !defined(AFX_UnitCell_H_8P4B21B8_C13C_5648_BF23_444095086234__INCLUDED_)
#define AFX_UnitCell_H_8P4B21B8_C13C_5648_BF23_444095086234__INCLUDED_


#include "SimDef.h"
#include "Vec3D.h"
class bead;
class UnitCell
{
public:
    

	UnitCell(int id, int i, int j, int k);
	 ~UnitCell();



        inline const std::vector <int>   GetIndex()    const    {return m_Index;}
        inline const int GetID()                const  		{return m_ID;}
        inline std::vector <bead *> GetBeadList()           {return m_BeadList;}
	inline Vec3D	GetCNTSidemin()              		{return m_CNTRmin;}
	inline Vec3D	GetCNTSidemax()              		{return m_CNTRmax;}
        inline std::vector <UnitCell *> GetVNeighbourUnitCell()   {return m_VNeighbourUnitCell;}




public:
    

  void AddtoVertexList(bead * z);
  void UpdateCNTSidemin(Vec3D);
  void UpdateCNTSidemax(Vec3D);
  void AddtoNeighbourUnitCell(UnitCell* z);
  void RemoveFromVertexList(bead * z);
    UnitCell * ReturnNewCell(int,int,int);

private:



private:


  int m_ID;
  std::vector <bead *> m_BeadList;
  std::vector <UnitCell *> m_VNeighbourUnitCell;
  std::vector <int>  m_Index;
  Vec3D		m_CNTRmin;
  Vec3D		m_CNTRmax;



    





};


#endif
