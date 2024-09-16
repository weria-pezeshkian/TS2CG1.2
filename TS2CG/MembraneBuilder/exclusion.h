#if !defined(AFX_exclusion_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_exclusion_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

class exclusion
{
public:
    
	exclusion(int id,int vid, double R);
	 ~exclusion();

	    inline const int GetID()            const  {return m_ID;}
        inline double GetRadius()               {return m_R;}
        inline const int GetPointID()       const  {return m_PointID;}

private:
    int m_PointID;
    int m_ID;
    double m_R;

};


#endif
