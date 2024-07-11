#if !defined(AFX_Wall_H_894B21B8_C13C_5648_BF23_124775086234__INCLUDED_)
#define AFX_Wall_H_894B21B8_C13C_5648_BF23_124775086234__INCLUDED_


#include "Def.h"
#include "Vec3D.h"
#include "UnitCell.h"
#include "point.h"
#include "bead.h"

class Wall
{
public:
    
	 Wall();
	 ~Wall();


    inline std::vector<bead> GetWallBead()                const  {return m_AllWallBeads;}
    inline std::vector<point*> GetWallPoint()                const  {return m_AllWallPoints;}
    inline bool GetState()                const  {return m_State;}
    inline double GetH()                const  {return m_H;}




public:
  void UpdateBox(Vec3D *x);
  void UpdateState(bool x);
  void UpdateUniform(bool x);
  void UpdateH(double x);
  void UpdateDen(double x);
  void UpdateBeadName(std::string x);
  void PrintWallState();
  void CreateWall(std::vector<point*>  p1, std::vector<point*>  p2);
  void Wall_itp();
    void UpdateCellSize(double x);

public:


private:
    bool m_Uniform;
    std::vector<bead> m_AllWallBeads;
    std::vector<point*> m_AllWallPoints;

  bool m_State;
  double m_Density;
  double m_H;
  std::string m_BeadName;
  Vec3D *m_pBox;
    double m_CellSize;
    std::vector<bead> MakeUniformBeads(std::vector<point*> &p);

    

};


#endif
