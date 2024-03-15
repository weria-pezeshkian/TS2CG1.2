#if !defined(AFX_Job_H_9F4A21B8_C13C_1223_BF23_124095086234__INCLUDED_)
#define AFX_Job_H_9F4A21B8_C13C_1223_BF23_124095086234__INCLUDED_

#include "SimDef.h"
class Job
{
public:
    
	Job(std::vector <std::string> argument);
	 ~Job();

private:
    std::string  TrimNameFromPath(std::string &); // small function to trim the path address from the name of the binary, we want to enforce fixed name, incase the binary was called in another software.
   std::vector <std::string> m_Argument;

};

#endif
