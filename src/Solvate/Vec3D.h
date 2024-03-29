#ifndef AFX_Vec3D_H_8F4421B8_D12D_11D3_CF24_124095086555__INCLUDED_
#define AFX_Vec3D_H_8F4421B8_D12D_11D3_CF24_124095086555__INCLUDED_

#include <vector>

class Vec3D {
public:
    Vec3D(double x, double y, double z);
    Vec3D();
    ~Vec3D();

public:
    double norm();
    double at(int n);
    void put(int n, double);

    double& operator()(const int n);
    Vec3D operator + (Vec3D);
    Vec3D operator - (Vec3D);
    Vec3D operator * (Vec3D);
    Vec3D operator * (double);

    void operator = (Vec3D);

    double dot(Vec3D, Vec3D);

private:
    double m_X;
    double m_Y;
    double m_Z;
    double m_Value;
};

#endif
