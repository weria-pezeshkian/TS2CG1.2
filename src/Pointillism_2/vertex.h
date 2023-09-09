#if !defined(AFX_vertex_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_vertex_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_


#include "SimDef.h"
#include "CNTCell.h"
#include "Vec3D.h"
#include "Tensor2.h"
#include "inclusion.h"

class links;
class triangle;
class vertex
{
public:
    
	vertex(int id, double x, double y, double z);
	vertex(int id);
	 ~vertex();


	inline const int GetVID()                const  {return m_ID;}

        // inline functions to access the privte veriables in the class
        inline double GetVXPos()                            {return m_X;}
        inline double GetVYPos()                            {return m_Y;}
        inline double GetVZPos()                            {return m_Z;}
        inline double GetArea()                             {return m_Area;}
        inline Tensor2  GetL2GTransferMatrix()              {return m_T_Local_2_Global;}
        inline Tensor2  GetG2LTransferMatrix()              {return m_T_Global_2_Local;}
        inline Vec3D GetNormalVector()                      {return m_Normal;}
        inline std::vector <double> GetCurvature()          {return m_Curvature;}
        inline std::vector <links *> GetVLinkList()             {return m_VLinkList;}
        inline std::vector <triangle *> GetVTraingleList()         {return m_VTraingleList;}
        inline std::vector <vertex *> GetVNeighbourVertex()     {return m_VNeighbourVertex;}
        inline CNTCell * GetVCNTCell()                      {return m_CNTCell;}
        inline inclusion* GetInclusion()                    {return m_pInclusion;}
        inline bool VertexOwnInclusion()                    {return m_OwnInclusion;}
        inline Vec3D *GetBox()                               {return m_pBox;}
        inline int GetDomainID()                               {return m_DomainID;}
        inline bool  GetIsFullDomain()                               {return m_IsFullDomain;}
        inline int GetClusterID()                      {return m_NetworkID;}



public:
    
  void UpdateVXPos(double x);    // a function for changing the x position of a vertex 
  void UpdateVYPos(double y);
  void UpdateVZPos(double z);
  void NOPBCUpdatePos(Vec3D z);
  void UpdateIsFullDomain(bool z);
  void UpdateVCNTCell(CNTCell * z);
  void AddtoLinkList(links* z);
  void AddtoTraingleList(triangle * z);
  void AddtoNeighbourVertex(vertex* z);
  void RemoveFromLinkList(links* z);
  void RemoveFromTraingleList(triangle * z);
  void RemoveFromNeighbourVertex(vertex* z);
  void UpdateBox(Vec3D *z);
  void UpdateCurvature(double,double); // vis
  void UpdateDomainID(int );                    // update domain ID, only in version 1.1 and above
  void UpdateNormal_Area(Vec3D,double); // vis
  void UpdateL2GTransferMatrix(Tensor2 v);
  void UpdateG2LTransferMatrix(Tensor2 v);
  void UpdateInclusion(inclusion* );
  void UpdateOwnInclusion(bool );
    void UpdateClusterID(int z);


public:
    bool CheckCNT();
    void ReadVertexFromFile(std::ifstream *inputfile,std::vector <vertex *> pv, std::vector <links *> pL, std::vector <triangle *> pT, std::vector <inclusion *> pI);
    void WriteVertexToFile(std::ofstream *inputfile);

private:

  int m_ID;
  double m_X;
  double m_Y;
  double m_Z;
    std::vector <triangle *> m_VTraingleList;
    std::vector <links *> m_VLinkList;
    std::vector <vertex *> m_VNeighbourVertex;
    inclusion *m_pInclusion;
    bool m_OwnInclusion;
    double m_Area;
    CNTCell * m_CNTCell;
    int m_DomainID;
    bool m_IsFullDomain;
    int m_NetworkID;            // which molecule or cluster the vertex belong too we will change this

////// Regarding the vertex curvature
private:


 Vec3D m_Normal;
 Vec3D *m_pBox;
 std::vector<double> m_Curvature;
    Tensor2  m_T_Local_2_Global;         //  Local to global transformation matrix
    Tensor2  m_T_Global_2_Local; 
};


#endif
