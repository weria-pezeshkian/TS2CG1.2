#include "Job.h"
#include "Nfunction.h"
#include "Edit_configuration.h"

Job::Job(std::vector <std::string> argument)
{
    /*
     *  Copyright Weria Pezeshkian (weria.pezeshkian@gmail.com), 2020.
     * This is a class to call different functions: at the moment, only one exist
     */

std::string exacutable=argument.at(0);
int  n=exacutable.size();
char L1 = exacutable.at(n-1);
char L2 = exacutable.at(n-2);
char L3 = exacutable.at(n-3);
bool condition = true;
	if (n>3)
	{
        char L4 = exacutable.at(n-4);
		if (L4!='/')
		{
			condition = false;
			std::cout<<argument.at(0)<<" 1. as a job name is wrong and does not exist \n";
		}

	}
    if(condition == true)
    {
        if( L3 == 'P' && L2 == 'L' && L1 == 'M')
        {

            Edit_configuration B(argument);
        }
        else if(n<3)
        {
            condition  = false;
            std::cout<<argument.at(0)<<"2. unrecognized executable binary name  \n";
        }

    }


}
Job::~Job()
{
    
}

