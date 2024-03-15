#include "Argument.h"
#include "help.h"
#include "Nfunction.h"

Argument::Argument(std::vector<std::string> argument) {
    m_Argument = argument;
    Nfunction f;

    std::string Arg1;
    m_In_GroFileName = "Input.gro";
    m_Out_GroFileName = "Output.gro";
    m_Tem_GroFileName = "W.gro";
    m_NegName = "CL";
    m_PosName = "NA";
    m_Seed = 9474;
    m_RCutOff = 0.4;
    m_DB = 0.05;
    m_UCELLSize = 2;
    m_Ion.push_back(0);
    m_Ion.push_back(0);

    // Iterate through the command-line options to assign values to the variables
    for (long i = 1; i < m_Argument.size(); i += 2) {
        Arg1 = m_Argument.at(i);
        if (Arg1 == "-in") {
            m_In_GroFileName = m_Argument.at(i + 1);
        } else if (Arg1 == "-h") {
            // Help message should be called
            help helpmessage(m_Argument.at(0));
            m_ArgCon = 0;
            exit(0);
        } else if (Arg1 == "-db") {
            m_DB = f.String_to_Double(m_Argument.at(i + 1));
        } else if (Arg1 == "-nname") {
            m_NegName = m_Argument.at(i + 1);
        } else if (Arg1 == "-pname") {
            m_PosName = m_Argument.at(i + 1);
        } else if (Arg1 == "-o") {
            m_Out_GroFileName = m_Argument.at(i + 1);
        } else if (Arg1 == "-ion") {
            int p = f.String_to_Double(m_Argument.at(i + 1));
            int n = f.String_to_Double(m_Argument.at(i + 2));
            m_Ion[0] = p;
            m_Ion[1] = n;
            i++;
        } else if (Arg1 == "-tem") {
            m_Tem_GroFileName = m_Argument.at(i + 1);
        } else if (Arg1 == "-seed") {
            m_Seed = f.String_to_Double(m_Argument.at(i + 1));
        } else if (Arg1 == "-unsize") {
            m_UCELLSize = f.String_to_Double(m_Argument.at(i + 1));
        } else if (Arg1 == "-Rcutoff") {
            m_RCutOff = f.String_to_Double(m_Argument.at(i + 1));
        } else {
            std::cout << "---> error: Wrong command: " << Arg1;
            std::cout << "\n" << "For more information and tips execute SOL -h" << "\n";
            m_ArgCon = 0;
            exit(0);
        }
    }
}
Argument::~Argument() {

}
