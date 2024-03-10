
#include "Argument.h"
#include "help.h"
#include "Nfunction.h"

Argument::Argument(std::vector <std::string> argument)
{
    m_Argument=argument;

    Nfunction f;
    
    

    
    
    m_Health             = true;
    m_TrjFile            = "MD.xtc";
    m_IniStep 			 = 0;
    m_FinStep 			 = 10000000;
    m_Seed                       = 7672;
    m_GeneralOutputFilename	 = "output";
	m_Jump = 1;
	m_AnalysisType= "Density_1D_2";
    m_GroFileName   =   "input.gro";
    m_SoftWareVersion = SoftWareVersion;
    m_IndexFileName = "Index.ndx";
    m_IndexFileName2 = "Index2.ndx";
    m_GrayAngle = 0;
    m_Cutoff = 1.1;
    m_MoreData = "no";
    m_NoBins = 100;
    m_TimeAverageFrames = 10;
    m_ArgCon=1;
    std::string Arg1;
    m_pdbFileName = "PDB.pdb";
    m_excutablename = m_Argument.at(0);
    m_binL = 0.1;
    m_TargetGroupName = "POPC";
    m_ReferenceGroupName = "PO4";
    
    if (m_Argument.size()>300)
    {
        std::cout << "Error: to many argument ";
        std::cout<<"\n"<<"For more information and tips execute MCL -h"<<"\n";
       m_ArgCon=0;
        m_Health = false;
        
    }
    else if (m_Argument.size() == 2)
    {
        
        Arg1 = m_Argument.at(1);
        if (Arg1 == "-h" )
        {
        // help message should be made
            help helpmessage(m_SoftWareVersion,m_excutablename);
            m_ArgCon=0;
            m_Health = false;
        }

    }
    else
    {

        for (long i=1;i<m_Argument.size();i=i+2)
        {
            Arg1 = m_Argument.at(i);
            if(Arg1=="-xtc")
            {
                m_TrjFile = m_Argument.at(i+1);
            }
            else if(Arg1=="-b")
            {
                m_IniStep = f.String_to_Int(m_Argument.at(i+1));
            }
            else if(Arg1=="-e")
            {
                m_FinStep = f.String_to_Int(m_Argument.at(i+1));
            }
            else if(Arg1=="-ndx")
            {
                m_IndexFileName = m_Argument.at(i+1);
            }
            else if(Arg1=="-tg")
            {
                m_TargetGroupName = m_Argument.at(i+1);
            }
            else if(Arg1=="-rg")
            {
                m_ReferenceGroupName = m_Argument.at(i+1);
            }
            else if(Arg1=="-gro")
            {
                m_GroFileName = m_Argument.at(i+1);
            }
            else if(Arg1=="-r")
            {
                m_AnalysisType = m_Argument.at(i+1);
            }
            else if(Arg1=="-j")
            {
                m_Jump = f.String_to_Int(m_Argument.at(i+1));
            }
            else if(Arg1=="-cutoff")
            {
                m_Cutoff = f.String_to_Double(m_Argument.at(i+1));
            }
            else if(Arg1=="-bl")
            {
                m_binL = f.String_to_Double(m_Argument.at(i+1));
            }
            else if(Arg1=="-more")
            {
                m_MoreData = m_Argument.at(i+1);
            }
            else if(Arg1=="-gangle")
            {
                m_GrayAngle = f.String_to_Double(m_Argument.at(i+1));
            }
            else if(Arg1=="-bins")
            {
                m_NoBins = f.String_to_Int(m_Argument.at(i+1));
            }
            else if(Arg1=="-tf")
            {
                m_TimeAverageFrames = f.String_to_Int(m_Argument.at(i+1));
            }
            else if(Arg1=="-seed")
            {
                m_Seed = f.String_to_Int(m_Argument.at(i+1));
            }
            else if(Arg1=="-defout")
            {
                m_GeneralOutputFilename = m_Argument.at(i+1);
            }
            else if(Arg1=="-ndx2")
            {
                m_IndexFileName2 = m_Argument.at(i+1);
            }
            else if(Arg1=="-pdb")
            {
                m_pdbFileName = m_Argument.at(i+1);
            }
            else
            {
                std::cout << "Error: wrong command :"<<Arg1;
                std::cout<<"\n"<<"For more information and tips execute ./DMC -h"<<"\n";
                m_ArgCon=0;
                m_Health = false;
            }
        }
        
    }


    
}






Argument::~Argument()
{
   
}
