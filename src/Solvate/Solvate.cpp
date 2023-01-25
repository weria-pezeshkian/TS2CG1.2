


#include "Argument.h"
#include "Solvate.h"
#include "Nfunction.h"
#include "GroFile.h"
#include "GenerateUnitCells.h"
#include "Solvate.h"

Solvate::Solvate(Argument *pArg)
{

    std::string ingrofilename  = pArg->GetIn_GroFileName();
    std::string outgrofilename = pArg->GetOut_GroFileName();
    std::string temfilename = pArg->GetTemplateFileName();
    double cutoff = pArg->GetRCutOff();
    int seed = pArg->GetSeed();
    std::vector<int> ion = pArg->GetIon();
    double db = pArg->GetDB();
    double usize  = pArg->GetUCELLSize();
    
    Nfunction f;
    bool health = true;
    if(f.FileExist(ingrofilename)==false)
    {
        std::cout<<"file "<<ingrofilename<< " do not exist. Error: at this moment we do not create empty box \n";
        health = false;
    }
    if(f.FileExist(temfilename)==false)
    {
        std::cout<<"file "<<temfilename<< " do not exist. Error: template file should be given  \n";
        health = false;
    }
    if(health==true)
    {
        GroFile InGro = GroFile(ingrofilename);
        
        std::vector<bead*> Sysbead = InGro.GetpAllBeads();

        Vec3D *FBox = InGro.GetBox();
        GenerateUnitCells UCELL(Sysbead,pArg,FBox, cutoff, usize);
        
        /// template water
        GroFile TemGro = GroFile(temfilename);
        std::vector<bead*> Wbead = TemGro.GetpAllBeads();
        Vec3D *WBox = TemGro.GetBox();
        
        int nBox_X = int((*FBox)(0)/(*WBox)(0))+1;
        int nBox_Y = int((*FBox)(1)/(*WBox)(1))+1;
        int nBox_Z = int((*FBox)(2)/(*WBox)(2))+1;
        std::vector<bead> FullWaterBead;

        for (int i=0;i<nBox_X;i++)
        {
            
        for (int j=0;j<nBox_Y;j++)
        {
        for (int k=0;k<nBox_Z;k++)
        {
            for (std::vector<bead *>::iterator it = Wbead.begin() ; it != Wbead.end(); ++it)
            {
                double x=(*it)->GetXPos()+((*WBox)(0))*double(i)+db;
                double y=(*it)->GetYPos()+((*WBox)(1))*double(j)+db;
                double z=(*it)->GetZPos()+((*WBox)(2))*double(k)+db;

                Vec3D Pos(x,y,z);
                if(UCELL.anythingaround (Pos)!=true)
                {

                if(x>0 && y>0 && z>0 && x<(*FBox)(0) && y<(*FBox)(1) && z<(*FBox)(2))
                {

                        bead TB = *(*it);
                
                        TB.UpdatePos(FBox,x,y,z);
                        FullWaterBead.push_back(TB);
                }
                }
            }
        }}}///  End full

        if(ion.at(0)!=0 || ion.at(1)!=0)
        {
            std::vector<bead> WaterBead;
            std::vector<bead> PIonBead;
            std::vector<bead> NIonBead;

            int tion = ion.at(0) + ion.at(1);
            
            int i=0;
            int noW = FullWaterBead.size();
            int P = noW/tion;
            int genion = 0;
            int genpion = 0;
            int gennion = 0;
            for (std::vector<bead>::iterator it = FullWaterBead.begin() ; it != FullWaterBead.end(); ++it)
            {
                if((i+1)%P==0 && genion<tion)
                {
                    double state = double(genpion-gennion)/(genion+1);
                    double target = double (ion.at(0)-ion.at(1))/tion;
                    if(state<target)
                    {
                        it->UpdateBeadName("Na");
                        it->UpdateResName("ION");
                        PIonBead.push_back(*it);

                        genpion++;

                    }
                    else if(gennion<ion.at(1))
                    {
                        it->UpdateBeadName("CL");
                        it->UpdateResName("ION");
                        NIonBead.push_back(*it);

                        gennion++;

                    }
                    else
                    {
                        it->UpdateBeadName("Na");
                        it->UpdateResName("ION");
                        PIonBead.push_back(*it);
                        genpion++;
                    }
                    
                    
                    genion++;
                }
                else
                {
                    WaterBead.push_back(*it);

                }
                
                i++;
            }
            
            std::cout<<" No of generated ions + "<<PIonBead.size()<<": - "<<NIonBead.size()<<"\n";
            std::cout<<" Requested + "<<ion.at(0)<<": - "<<ion.at(1)<<"\n";

            std::ofstream info;
            info.open("info.txt");
            info<<"W    "<<WaterBead.size()<<"\n";
            info<<"NA    "<<PIonBead.size()<<"\n";
            info<<"CL    "<<NIonBead.size()<<"\n";
             std::vector<bead> PreBeads = InGro.GetAllBeads();
            for (std::vector<bead>::iterator it = WaterBead.begin() ; it != WaterBead.end(); ++it)
                PreBeads.push_back(*it);
            
            for (std::vector<bead>::iterator it = PIonBead.begin() ; it != PIonBead.end(); ++it)
                PreBeads.push_back(*it);
            
            for (std::vector<bead>::iterator it = NIonBead.begin() ; it != NIonBead.end(); ++it)
                PreBeads.push_back(*it);
            
            TemGro.RenewBeads(PreBeads);
            TemGro.UpdateBox(*FBox);
            TemGro.WriteGroFile(outgrofilename);
            
            
            
        }
        else
        {
        
            std::vector<bead> PreBeads = InGro.GetAllBeads();
            for (std::vector<bead>::iterator it = FullWaterBead.begin() ; it != FullWaterBead.end(); ++it)
                PreBeads.push_back(*it);

            TemGro.RenewBeads(PreBeads);
            TemGro.UpdateBox(*FBox);
            TemGro.WriteGroFile(outgrofilename);
            
            std::ofstream info;
            info.open("info.txt");
            info<<"W    "<<FullWaterBead.size()<<"\n";

        }
        
        
    }
    



}
Solvate::~Solvate()
{
    
}

