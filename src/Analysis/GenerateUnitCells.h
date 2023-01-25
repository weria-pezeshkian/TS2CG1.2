#if !defined(AFX_GenerateUnitCells_H_8P4B21B8_C13C_5648_BF23_444095086239__INCLUDED_)
#define AFX_GenerateUnitCells_H_8P4B21B8_C13C_5648_BF23_444095086239__INCLUDED_


#include "SimDef.h"
#include "UnitCell.h"
#include "Argument.h"
#include "bead.h"
#include "Vec3D.h"
class GenerateUnitCells
{
public:
    

	GenerateUnitCells(std::vector< bead* > bead,Argument *pArgu,Vec3D *pBox);
	~GenerateUnitCells();




        inline std::vector <UnitCell *> GetAllCNTCells()           {return m_pAllCNTCells;}
    inline std::vector <double> GetCNTCellSize()        {return m_CNTCellSize;}
    inline std::vector <int> GetCNTCellNo()        {return m_CNTCellNo;}




public:
    

//  void AddtoVertexList(vertex * z);

    void Generate();

private:
    Argument *m_pArgu;
std::vector <UnitCell *> m_pAllCNTCells;
std::vector <UnitCell > m_AllCNTCells;
int IDFromIndex(int,int,int);
int IndexFromID(int,int *,int *,int *);
double m_CNTSize;
private:
    std::vector< bead* > m_pAllBead;

int m_Nx;
int m_Ny;
int m_Nz;
    Vec3D  *m_pBox;
    std::vector <double> m_CNTCellSize;
    std::vector <int> m_CNTCellNo;




    





};


#endif
