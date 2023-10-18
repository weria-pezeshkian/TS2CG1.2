  #if !defined(AFX_GenDomains_CPP_334A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_GenDomains_CPP_334A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "GenDomains.h"
#include "Def.h"

/*
 

 
 
 
 */

GenDomains::GenDomains(std::string strfilename, std::vector<point*>  point1,std::vector<point*>  point2,bool renorm)
{
    m_Health = true;
    m_Point1 = point1;
    m_Point2 = point2;
    
    //************************ Read str file *************/
    
    Nfunction f;
    std::vector <int> DomainTypes;
    std::ifstream strfile;
    strfile.open(strfilename.c_str());
    std::string str,lname;
    double ratio1,ratio2;
    double ap;
    bool flag =false;
    int FileLine = 0;
    while (true)
    {
        std::getline (strfile,str);
        FileLine++;
        if(strfile.eof())
        break;
        

        std::vector<std::string> Line = f.split(str);

    if(Line.size()!=0 && (Line.at(0)).at(0)!=';')
    {
        if((Line.at(0)).at(0)=='[' )
        {
            str = f.trim(str);
            str.erase(std::remove(str.begin(), str.end(), '['), str.end());
            str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
            str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
            if(str=="LipidsList")
            {
                flag = true;
            }
            else if(str=="ProteinList")
            {
                flag = false;
            }
            else if(str=="ShapeData")
            {
                flag = false;
            }
            else
            {
                std::cout<<" Error: unknown command in line "<<FileLine<<" <"<<str<<"> of the "<<strfilename<< " file \n";
                m_Health=false;
            }
        }
        else if(flag==true )
        {
            if(Line.at(0)=="Domain" && Line.size()>1)
            {

                int domainid = f.String_to_Int(Line.at(1));
                //std::getline (strfile,str);
                std::vector<point*> P1,P2;
                DomainTypes.push_back(domainid);
                for ( std::vector<point*>::iterator it = m_Point1.begin(); it != m_Point1.end(); it++ )
                {
                    int id = (*it)->GetDomainID();
                    if(id==domainid)
                    P1.push_back((*it));
                }
                for ( std::vector<point*>::iterator it = m_Point2.begin(); it != m_Point2.end(); it++ )
                {
                    int id = (*it)->GetDomainID();
                    if(id==domainid)
                    P2.push_back((*it));
                }
                Domain Do1 (domainid,P1);
                Domain Do2 (domainid,P2);

                    while (true)
                    {
                            FileLine++;
                            std::getline (strfile,str);
                            std::vector<std::string> LL = f.split(str);
                            if(LL.size()!=0 && (LL.at(0)=="END" || LL.at(0)=="End"))
                                break;
                        
                            if(LL.size()>=4 && (LL.at(0)).at(0)!=';')
                            {
                                Do1.AddADomainLipid(LL.at(0), f.String_to_Double(LL.at(3)), f.String_to_Double(LL.at(1)));
                                Do2.AddADomainLipid(LL.at(0), f.String_to_Double(LL.at(3)), f.String_to_Double(LL.at(2)));
                                if(f.String_to_Double(LL.at(3))<0 || f.String_to_Double(LL.at(2))<0 ||f.String_to_Double(LL.at(1))<0)
                                {
                                    std::cout<<" Error: line "<<FileLine<<" of the file "<<strfilename<<"; contains negative values  \n";
                                    m_Health=false;
                                }
                            }

                    }

                m_AllDomains.push_back(Do1);
                m_AllDomains.push_back(Do2);


            }
            else
            {
                std::cout<<" Error 232----GENDOMAIN:   "<<flag<<"  "<<str<<"\n";

            }

        }
    }
    }

    strfile.close();

    for ( std::vector<point*>::iterator it = m_Point2.begin(); it != m_Point2.end(); it++ )
    {
        int id = (*it)->GetDomainID();
        bool found = false;
        for ( std::vector<int>::iterator it2 = DomainTypes.begin(); it2 != DomainTypes.end(); it2++ )
        {
            if(id==(*it2))
            found = true;
        }
        if(found == false)
        {
            std::cout<<" Error: there are points with domain id of "<<id<<" while this id is not defined in the str file \n";
            std::exit(0);

        }

        
    }
    /// CHECK if all the points domain are defined in the file
    
    for ( std::vector<Domain>::iterator it = m_AllDomains.begin(); it != m_AllDomains.end(); it++ )
    {
        m_pAllDomains.push_back(&(*it));
    }
    for ( std::vector<Domain*>::iterator it = m_pAllDomains.begin(); it != m_pAllDomains.end(); it++ )
        (*it)->Configure(renorm);


}
GenDomains::~GenDomains()
{
    
}





#endif



