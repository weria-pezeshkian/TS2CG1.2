


#include "Argument.h"
#include "Job.h"
#include "Nfunction.h"
#include "Curvature_Corr.h"
#include "TwoLipid_Distribution.h"
#include "Curvature_Hist.h"
#include "Sin_1D.h"
#include "Density_1D.h"
#include "Density_1D_2.h"


Job::Job(std::vector <std::string> argument)
{



Argument a(argument);
    std::string jobtype = a.GetJobType();



if(a.GetHealth()==true )
{

	if( jobtype == "corr"  )
	{

		Curvature_Corr B(&a);

	}
    if( jobtype == "2com"  )
    {
        
        TwoLipid_Distribution B(&a);
        
    }
    if( jobtype == "hist"  )
    {
        
        Curvature_Hist B(&a);
        
    }
    if( jobtype == "sin1D"  )
    {
        
        Sin_1D B(&a);
        
    }
    if( jobtype == "Density_1D"  )
    {
        
        Density_1D B(&a);
        
    }
    if( jobtype == "Density_1D_2"  )
    {
        
        Density_1D_2 B(&a);
        
    }
	else
	{
		std::cout<<jobtype<<" has not been implemented yet \n";
	}
}
else
{
std::cout<<" Error in input files \n";
}

}
Job::~Job()
{
    
}

