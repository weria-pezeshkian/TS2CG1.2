#if !defined(AFX_bead_H_334B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_bead_H_334B21B8_C13C_5648_BF23_124095086234__INCLUDED_
/* A bead object encapsulates the following attributes:
1) ID 2) Name 3) Type 4) Residue name 5-7) Position coordinates 8) Residue ID 9) A pointer to a unit cell.
The unit cell pointer enables efficient neighbor bead queries. While it's possible to perform such queries without the unit cell pointer, having it as a member function enhances convenience. */
#include <string>
// Forward declaration to reduce unnecessary inclusion of headers
class Vec3D;
class UnitCell;

class bead {
public:
    // Constructors and destructor
    bead(int id, const std::string& name, const std::string& type, const std::string& resname, int resid, double x, double y, double z);
    bead(int id, const std::string& name, const std::string& type, const std::string& resname, int resid);
    ~bead();

    // Accessor functions
    inline const int GetID() const { return m_ID; }
    inline double GetXPos() const { return m_X; }
    inline double GetYPos() const { return m_Y; }
    inline double GetZPos() const { return m_Z; }
    inline const std::string& GetBeadType() const { return m_BeadType; }
    inline const std::string& GetResName() const { return m_ResName; }
    inline int GetResid() const { return m_Resid; }
    inline const std::string& GetBeadName() const { return m_BeadName; }
    inline UnitCell* GetBeadUnitCell() const { return m_BeadUnitCell; }

    // Update functions
    void UpdateXPos(double x);
    void UpdateYPos(double y);
    void UpdateZPos(double z);
    void UpdateBox(Vec3D* z);
    void UpdateBeadName(const std::string& name);
    void UpdateResName(const std::string& name);
    void UpdatePos(Vec3D* B, double x, double y, double z);
    void UpdatePos(double x, double y, double z);
    void UpdateBeadUnitCell(UnitCell* z);

private:
    bool m_hasMol;
    double m_X;
    double m_Y;
    double m_Z;
    std::string m_BeadName;
    std::string m_BeadType;
    std::string m_ResName;
    int m_Resid;
    int m_ID;
    Vec3D* m_pBox;
    UnitCell* m_BeadUnitCell;
};

#endif
