#if !defined(AFX_Job_H_9F4A21B8_C13C_1223_BF23_124095086234__INCLUDED_)
#define AFX_Job_H_9F4A21B8_C13C_1223_BF23_124095086234__INCLUDED_

#include "Def.h"
class Job
{
public:
    
	Job(std::vector <std::string> argument);
	 ~Job();
    


public:
    
       // inline const std::string GetFileName()    const {return m_FileName;}




private:
    std::string SubstringFromRight(const std::string& input);

   std::vector <std::string> m_Argument;

};

#endif
