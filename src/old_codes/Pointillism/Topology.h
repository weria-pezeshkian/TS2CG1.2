
#if !defined(AFX_Topology_H_7F4B21B8_D13C_9321_QF23_124095086234__INCLUDED_)
#define AFX_Topology_H_7F4B21B8_D13C_9321_QF23_124095086234__INCLUDED_
/*
*  Copyright Weria Pezeshkian (weria.pezeshkian@gmail.com), 2020.
* To read and write q files.
*/
#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "Vec3D.h"

class Topology
{
public:
    
    Topology(Vec3D *Box ,double * pminangle, double * pLmin , double * pLmax);
    ~Topology();
    
    inline std::vector<vertex*> GetVertex()                     {return m_pAllV;} // returns pointers of all vertices
    inline std::vector<triangle*> GetTriangle()                 {return m_pAllT;} // returns pointera of all triangles
    inline std::vector<links*> GetLinks()                       {return m_pLinks;} // returns pointers of all links (note: v1-v2 and v2-v1 are two different links)
    inline std::vector<links*> GetHalfLinks()                    {return m_pHL;} // half of the links
    inline std::vector<links*> GetMHalfLinks()                   {return m_pMHL;} // the other half
    inline const  bool GetTopologyHealth()                       {return m_TopologyHealth;} // to check if any error has ever happened within this class
    
public:
    
    void Write(std::vector<vertex* > ver, std::vector<triangle* > triangle,  std::string Filename); //Function to  write a *.q file
    void FastReadQFile(std::string input ); //Function to  read a *.q file

private:
    
    std::vector<vertex> m_Vertex; // Keeps the orginal objects of the vertices
    std::vector<triangle> m_Triangle;
    std::vector<links> m_Links;
private:
    
    std::vector<vertex*>      m_pAllV; // pointers of all vertices, to be shared with other classes
    std::vector<triangle*>    m_pAllT;
    std::vector<links*>       m_pLinks;
    std::vector<links*>       m_pHL;
    std::vector<links*>       m_pMHL;
    
    Vec3D m_Box;        // System box size
    Vec3D *m_pBox;      // a pointer to the system box size,
    
    
    double *m_pLmin2;   // not relevant for this code, but for dts simulation, it should not be larger than L
    double *m_pLmax2;   // not relevant for this code, but for dts simulation, it should not be larger than sqrt(3)L
    double *m_pminAngle; // not relevant for this code, but for dts simulation is
    bool m_TopologyHealth; //returns false whenever an error happened in this class
private:
    double LengthBetweenTwoVertex(vertex* v1, vertex* v2);   // distance between to vertices
    bool   CheckFaceAngle(links * l);  // calculate the angle between two neighboring triangle
    Vec3D   CalculateNormal(vertex*,vertex*,vertex*); // assumes v1,v2,v3 are three vertices of a trinagle and gives notmal vector
    
};


#endif
