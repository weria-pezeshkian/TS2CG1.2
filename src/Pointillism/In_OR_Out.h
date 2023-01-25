#if !defined(AFX_In_OR_Out__INCLUDED_)
#define AFX_In_OR_Out__INCLUDED_
/*
 *  Copyright Weria Pezeshkian (weria.pezeshkian@gmail.com), 2020.
 * This is a class to call different functions: at the moment, only one exist
 */
#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
/*
 this class is takes a TS file, ask for a point and checks if the point is inside the TS file.
 */
class In_OR_Out
{
public:
    
	In_OR_Out(std::string);
	 ~In_OR_Out();

private:

Vec3D *m_pBox;
    Vec3D m_Box;

    std::vector<vertex*>      m_pAllV;
    std::vector<vertex>       m_AllV;
    std::vector<triangle*>    m_pAllT;
    std::vector<links*>       m_pAllLinks;
    std::vector<links*>       m_pHalfLinks1;
    std::vector<links*>       m_pHalfLinks2;
    std::vector<inclusion*>   m_pInc;
    void UpdateGeometry( );  // updates curvature, area etc of each triangle, vertex etc

    void Initialize(std::string file);
    bool FileExist (const std::string& name);




};


#endif
