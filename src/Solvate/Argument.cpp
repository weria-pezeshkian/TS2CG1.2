
#include "Argument.h"
#include "help.h"
#include "Nfunction.h"

Argument::Argument(std::vector <std::string> argument)
{
    m_Argument=argument;

    Nfunction f;
    
    

    std::string Arg1;
    m_In_GroFileName   =  "Input.gro";
    m_Out_GroFileName  =  "Output.gro";
    m_Tem_GroFileName  = "W.gro";
    m_Seed             = 9474;
    m_RCutOff          = 0.4;
    m_Health             = true;
    m_DB                = 0.05;
    m_UCELLSize         = 2;
    m_Ion.push_back(0);
    m_Ion.push_back(0);

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
            help helpmessage(m_Argument.at(0));
            m_ArgCon=0;
            m_Health = false;
        }

    }
    else
    {

        for (long i=1;i<m_Argument.size();i=i+2)
        {
            Arg1 = m_Argument.at(i);
            if(Arg1=="-in")
            {
                m_In_GroFileName = m_Argument.at(i+1);
            }
            else if(Arg1=="-db")
            {
                m_DB = f.String_to_Double(m_Argument.at(i+1));
            }
            else if(Arg1=="-o")
            {
                m_Out_GroFileName = m_Argument.at(i+1);
            }
            else if(Arg1=="-ion")
            {
                int p =  f.String_to_Double(m_Argument.at(i+1));
                int n =  f.String_to_Double(m_Argument.at(i+2));
                m_Ion.at(0) = p;
                m_Ion.at(1) = n;
                i++;

            }
            else if(Arg1=="-tem")
            {
                m_Tem_GroFileName = m_Argument.at(i+1);
            }
            else if(Arg1=="-seed")
            {
                m_Seed = f.String_to_Double(m_Argument.at(i+1));
            }
            else if(Arg1=="-unsize")
            {
                m_UCELLSize = f.String_to_Double(m_Argument.at(i+1));
            }
            else if(Arg1=="-Rcutoff")
            {
                m_RCutOff = f.String_to_Double(m_Argument.at(i+1));
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
