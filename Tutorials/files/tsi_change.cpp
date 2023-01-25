#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <list>
#include <map>
#include <iomanip>
#include <valarray>
#include <stdio.h>
struct vertex {
    
    int id;
    double x;
    double y;
    double z;
    int domain;
    
} ;
struct triangle {
    
    int id;
    int v1;
    int v2;
    int v3;
    
} ;
struct inclusion {
   
    int id;
    int type;
    int v;
    double x;
    double y;
    
} ;
struct tsi {
   
    std::vector<vertex> ver;
    std::vector<triangle> tri;
    std::vector<inclusion> inc;
    double lx;
    double ly;
    double lz;
    
} ;
tsi readtsi(std::string);
std::vector<std::string> split(std::string str);
void writetsi(tsi TSI,std::string filename);
double S2D(std::string ConInt);
int main(int argc, char* argv[])
{
    std::vector <std::string> argument;
    std::string Temprory;
           for (long i=0;i<argc;i++)
           {
               Temprory.assign(argv[i]);
               argument.push_back(Temprory);
           }
if(argc<4)
{
std::cout<<"error: no tsi input and output-- two arguments are needed \n";
}
{
tsi TSI = readtsi(argument.at(1));
std::cout<<" we read the vertices \n";

    std::vector<vertex> ver= TSI.ver;
    std::vector<inclusion> inc =TSI.inc;
	double D = S2D(argument.at(3));
    for (int i=0;i<inc.size();i++)
	{
		vertex V = ver.at(((inc.at(i)).v));

		    for (int j=0;j<ver.size();j++)
			{
				double x2 = (V.x-(ver.at(j)).x)*(V.x-(ver.at(j)).x);
				double y2 = (V.y-(ver.at(j)).y)*(V.y-(ver.at(j)).y);
				double z2 = (V.z-(ver.at(j)).z)*(V.z-(ver.at(j)).z);
				double d2=x2+y2+z2;
				if(d2<D*D)
				(ver.at(j)).domain = (inc.at(i)).type;
			}

	}

TSI.ver = ver;
TSI.inc = inc;
std::cout<<" now we are writing the final tsi file \n";
writetsi(TSI,argument.at(2));
}

        








}
tsi readtsi(std::string filename)
{

    tsi TSI;
    std::string str;
    double x,y,z;
    int j,k,n,m,NO;
    std::ifstream input;
    input.open(filename.c_str());
    input>>str>>x>>str>>x>>y>>z;
    TSI.lx=x;
    TSI.ly=y;    
    TSI.lz=z;
    input>>str>>NO;
    getline(input,str);
    std::vector<vertex> ver;
    for (int i=0;i<NO;i++)
	{
		vertex V;
    		getline(input,str);
		std::vector<std::string> S = split(str);
		V.id  = S2D(S.at(0));
		V.x  = S2D(S.at(1));
		V.y  = S2D(S.at(2));
		V.z  = S2D(S.at(3));
		V.domain  = 0;
		if(S.size()>4)
		V.domain  = S2D(S.at(4));
		ver.push_back(V);

	}
		TSI.ver = ver;
    input>>str>>NO;
    std::vector<triangle> tri;
    for (int i=0;i<NO;i++)
	{
   		 input>>j>>k>>n>>m;
		triangle TRI;
		TRI.id  = i;
		TRI.v1  = k;
		TRI.v2  = n;
		TRI.v3  = m;
		tri.push_back(TRI);

	}
		TSI.tri = tri;
    input>>str>>NO;
    std::vector<inclusion> inc;
    for (int i=0;i<NO;i++)
	{
   		input>>j>>k>>n>>x>>y;
		inclusion INC;
		INC.id  = i;
		INC.type  = k;
		INC.v  = n;
		INC.x  = x;
		INC.y  = y;
		inc.push_back(INC);
	}
		TSI.inc = inc;





return TSI;
}
void writetsi(tsi TSI,std::string filename)
{

    FILE * output;
    output = fopen(filename.c_str(), "w");
    const char* versiontsi="version 1.1";
    fprintf(output,"%s\n",versiontsi);
    const char* box="box";
    fprintf(output,"%s%18.10lf%18.10lf%18.10lf\n",box,TSI.lx,TSI.ly,TSI.lz);
    const char* ver="vertex";
    int size=(TSI.ver).size();
    fprintf(output,"%s%20d\n",ver,size);

   for (int i = 0;i<size;i++)
        fprintf(output,"%d%18.10lf%18.10lf%18.10lf%10d\n",i,((TSI.ver).at(i)).x,((TSI.ver).at(i)).y,((TSI.ver).at(i)).z,((TSI.ver).at(i)).domain);

    const char* tri="triangle";
    size = (TSI.tri).size();
    fprintf(output,"%s%20d\n",tri,size);


   for (int i = 0;i<size;i++)
        fprintf(output,"%10d%10d%10d%10d\n",i,((TSI.tri).at(i)).v1,((TSI.tri).at(i)).v2,((TSI.tri).at(i)).v3);

    const char* inc="inclusion";
    size = (TSI.inc).size();
    fprintf(output,"%s%20d\n",inc,size);

   for (int i = 0;i<size;i++)
        fprintf(output,"%10d%10d%10d%18.10lf%18.10lf\n",i,((TSI.inc).at(i)).type,((TSI.inc).at(i)).v,((TSI.inc).at(i)).x,((TSI.inc).at(i)).y);




}
std::vector<std::string> split(std::string str)
{
    
    std::vector<std::string> Line;
    
    std::string word;
    bool flag=false;
    for (int i=0;i<str.size();i++)
    {
        if (str.at(i) ==' '||str.at(i) =='\f' || str.at(i) =='\n' || str.at(i) =='\r' || str.at(i) =='\t' || str.at(i) =='\v' )
        {
            if(flag == true)
            {
                Line.push_back(word);
                word.clear();
                flag=false;
            }
        }
        else
        {
            flag = true;
            word.push_back(str.at(i));
        }
    }
    if(flag == true)
    {
        Line.push_back(word);
        word.clear();
        flag=false;
    }
    return Line;
}
double S2D(std::string ConInt)
{
 
     double i =atof( ConInt.c_str() );
    return i;
}
