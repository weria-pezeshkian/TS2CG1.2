import numpy as np
import sys
import traceback
from inspect import cleandoc as cd
import os,shutil
import argparse


class PointUpdaterClass():
    """
    The point class is the class equivalent to the point folder created by TS2CG. The accepted parameters are discussed in the corresponding __init__.
    In short, each of the four files receives its own subclass. Additionally, some manipulating and saving routines are provided.
    """

    def _get_box(self,file):
        with open(file,"r",encoding="UTF8") as f:
            lines=f.readlines()[:4]
            for line in lines:
                if "Box" in line:
                    self.box=line
                    break


    def read_BM(self,file,coords=True):
        """
        read_BM reads the files in the point folder to return it to the class.
        :param file: the path to the file in the point folder
        :param coords: (default=True) For the Inclusions and Exclusions files (not the inner and outer layer), coords is set to False

        :return: the data contained in the file in a transposed manner to be used in the appropriate subclasses.
        """
        if coords:
            self._get_box(file)
            try:
                data=np.loadtxt(file,skiprows=3)
                return data.T
            except ValueError:
                data=np.loadtxt(file,skiprows=4)
                return data.T
        else:
            try:
                data=np.loadtxt(file,skiprows=2)
                return data.T
            except (ValueError,FileNotFoundError):
                return None

    class InnerOrOuter():
        """
        The InnerOrOuter class provides the framework to read InnerBM.dat and OuterBM.dat
        """

        def make_dict(self):
            names="id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2".split()
            for i,row in enumerate(self.raw):
                    self.get_data[names[i]]=row
            for name in ["id","domain_id"]:
                self.get_data[name]=np.asarray(self.get_data[name],dtype=int)

        def __init__(self,data=None,h_Name="",init=False):
            """
            :param data: The data of the appropriate dat file as has been read by read_BM()
            """
            
            self.h_Name=h_Name
            self.get_data={}
            if not init:
                if data is None:
                    print(traceback.format_exc())
                    sys.exit(cd("No data was passed to the Inner our Outer segment. Nothing can be done."))
                self.raw=data
                self.make_dict()
            else:
                self.get_data=dict.fromkeys(names,None)


    class Inclusions():
        """
        Represents the IncData.dat from the point folder.
        """
        def make_data(self,raw):
            names="id typeid pointid lx ly lz".split()
            self.get_data={}
            try:
                for i,row in enumerate(raw):
                    self.get_data[names[i]]=row
            except TypeError:
                pass

        def __init__(self,data):
            """
            :param data: The data of the appropriate dat file as has been read by read_BM()
            """
            self.raw=data
            try:
                self.NoInc=len(data)
            except TypeError:
                self.NoInc=0
            self.make_data(self.raw)




    class Exclusions():
        """
        Represents the ExcData.dat from the point folder.
        """
        def __init__(self,data):
            """
            :param data: The data of the appropriate dat file as has been read by read_BM()
            """
            names="id  pointid r".split()
            self.get_data={}
            self.raw=data
            try:
                for i,row in enumerate(self.raw):
                    self.get_data[names[i]]=row
            except TypeError:
                pass


    def __init__(self, path="point/",k=1,new=False):
        """
        Initliazes the point class to represent the point folder in the python framework. It creates the necessary
        subclasses for the whole folder named
            inner
            outer
            inclusions
            exclusions
        and initializes a place holder for the input_file which should have the form (for example):
            ; domain lipid percentage c0 density
            0 POPC .3 0.179 0.64
            1 POPD .4 0.613 0.64
            2 POPG .3 0.629 0.64
        depending on the alteration planned later.
        If necessarry k can be adjusted to correct the exponent in the assign_by_c0 formula.

        :param path: (default="point/)") gives the path to the point folder to be read out.
        """
        self.box="Box XXX XXX XXX"
        if path[-1]=="/":
            self.path=path
        else:
            self.path=path+"/"
        try:
            self.inner=self.InnerOrOuter(self.read_BM(self.path+"InnerBM.dat"),h_Name="inner")
            self.outer=self.InnerOrOuter(self.read_BM(self.path+"OuterBM.dat"),h_Name="outer")
            self.inclusions=self.Inclusions(self.read_BM(self.path+"IncData.dat",coords=False))
            self.exclusions=self.Exclusions(self.read_BM(self.path+"ExcData.dat",coords=False))
        except:
            print(traceback.format_exc())
            sys.exit(cd("Something went wrong initializing the point class. Please make sure that\n\
            * path is set to the point folder\n\
            * the point folder contains the files InnerBM.dat, OuterBM.dat, IncData.dat, and ExcData.dat"))
        self.input_file=None
        self.k=k
        self.protein_altered=False

    def update_domains(self,**kwargs):
        """
        update_domains is a helper function to manually update the domains in the point folder to place lipids. The function can be
        completely replaced by altering inner and outer directly from the point class. The method here merely provides some sanity checks
        and errors out, if it cannot be correct what was done.

        :param layer: (default=both) a parameter that chooses if the *upper*, *lower*, or *both* layers are altered.
        :param method: (default=None) a function can be passed to alter the domains
        :param assigned: (default=None) takes a numpy array with the same length as points to be assigned to domains and can
        be used to directly set the domains according to a list.

        """
        layer=kwargs.get("layer","both").lower()
        method=kwargs.get("method",None)
        assigned=kwargs.get("assign",None)
        domains={}
        if layer=="both":
            domains["Inner"]=self.inner.get_data["domain_id"]
            domains["Outer"]=self.outer.get_data["domain_id"]
        elif layer=="inner":
            domains["Inner"]=self.inner.get_data["domain_id"]
        elif layer=="outer":
            domains["Outer"]=self.outer.get_data["domain_id"]
        else:
            print(''.join(traceback.format_stack()))
            sys.exit(cd(f"update_domains allows either to change\n\
                * both\n\
                * inner\n\
                * outer\n\
                '{layer}' is not a valid choice."))


        if method is None and assigned is None:
            print("No change is passed to update_domains")
        elif assigned is not None:
            for key in domains:
                _length=len(domains[key])
                if len(assigned)==_length:
                    domains[key]=np.asarray(assigned)
                else:
                    print(''.join(traceback.format_stack()))
                    sys.exit(cd(f"update_domains only allows to replace the domains with a domain table\n\
                        of the same size.\n\
                        A table was supplied with size {len(assigned)}, but {_length} needs to fit."))
        elif method is not None:
            for key in domains:
                _length=len(domains[key])
                new_domain=method(**kwargs)
                if _length!=len(new_domain):
                    print(''.join(traceback.format_stack()))
                    sys.exit(cd(f"update_domains only allows to replace the domains, if the replacement\n\
                        has the same length. The supplied method did not return a correct length list.ņ\
                        It returned a list with size {len(new_domain)}, but it needs to be {_length}."))
                domains[key]=new_domain

        try:
            self.inner.get_data["domain_id"]=domains["Inner"]
        except KeyError:
            pass
        try:
            self.outer.get_data["domain_id"]=domains["Outer"]
        except KeyError:
            pass

    def assign_by_c0(self,input_file,location="both",unspecified_Number=False):
        """
        the function assign_by_c0 reads an input file of form (axample)
            ; domain lipid percentage c0 density
            0 POPC .3 0.179 0.64
            1 POPD .4 0.613 0.64
            2 POPG .3 0.629 0.64
        which would set domains 0, 1, 2 for the lipids randomly with their corresponding percentage such that the given c0 is closest to
        still unoccupied places.

        :param input_file: path to the input file
        :param location: (default="both") allows the alteration of just the *upper*,*lower*, or *both* monolayer.
        """
        with open(input_file,"r",encoding="UTF8") as f:
            lines=f.readlines()
            lines=[line for line in lines if line.strip()]
            lines=[line[:line.find(";")].split() for line in lines if line[0] != ";"]
            if location=="outer":
                N=len(self.outer.get_data["id"])
                locations=[self.outer]
            elif location=="inner":
                N=len(self.inner.get_data["id"])
                locations=[self.inner]
            elif location=="both":
                N=len(self.inner.get_data["id"])
                locations=[self.inner,self.outer]
            else:
                print(''.join(traceback.format_stack()))
                sys.exit(cd(f"assign_by_c0 allows either to consider\n\
                    * both\n\
                    * inner\n\
                    * outer\n\
                    '{location}' is not a valid choice."))


            if len(self.inner.get_data["id"])>=len(self.outer.get_data["id"]):
                randomizer=np.random.permutation(self.inner.get_data["id"])
            else:
                randomizer=np.random.permutation(self.outer.get_data["id"])

            for loc in locations:
                if loc.h_Name=="inner":
                    turn=-1
                else:
                    turn=1
                lipids={}
                for item in lines:
                    lipids[item[0]]=[round(float(item[2])*N,0)+1,float(item[3])]
                for index in randomizer:
                    domain=0
                    lipid_probabilities={}
                    Cs=np.asarray([loc.get_data["C1"][index],loc.get_data["C2"][index]])
                    for key in lipids:
                        Cs_input=np.asarray([lipids[key][1]])[0]
                        lipid_probabilities[key]=np.exp(-self.k*(turn*(Cs[0]+Cs[1])-Cs_input)**2)

                    normalizer=np.sum(np.asarray(list(lipid_probabilities.values())))
                    if normalizer==0:
                        for key in lipid_probabilities:
                            lipid_probabilities[key]=1/len(lipid_probabilities.keys())
                            normalizer=1


                    for key in lipid_probabilities:
                        lipid_probabilities[key]=lipid_probabilities[key]/normalizer
                    domain=np.random.choice(list(lipid_probabilities.keys()),1,p=list(lipid_probabilities.values()))[0]
                    deleter=None
                    for key in lipids:
                        if int(domain)==int(key) and not unspecified_Number:
                            lipids[key][0]=lipids[key][0]-1
                            if lipids[key][0]==0 and len(lipids.keys())>1:
                                deleter=key
                    if deleter is not None:
                        del lipids[deleter]
                        deleter=None
                    loc.get_data["domain_id"][index]=domain

    def assign_by_function(self,input_file,method,location="both",unspecified_Number=False):
        """
        the function assign_by_function assigns placement probabilities based on a method. The method only needs to be able
        to return a number per lipid.

        :param input_file: The input file contains the number per lipid type, etc. And might be a good place to define
        additional parameters for the method.
        :param method: a method that sorts the sorting probabilities
        :param location: (default="both") allows the alteration of just the *upper*,*lower*, or *both* monolayer.
        :param unspecified_Number: (default=False) If set, it will disregard
        
        """
        with open(input_file,"r",encoding="UTF8") as f:
            lines=f.readlines()
            lines=[line[:line.find(";")].split() for line in lines if line[0] != ";"]
            if location=="outer":
                N=len(self.outer.get_data["id"])
                locations=[self.outer]
            elif location=="inner":
                N=len(self.inner.get_data["id"])
                locations=[self.inner]
            elif location=="both":
                N=len(self.inner.get_data["id"])
                locations=[self.inner,self.outer]
            else:
                print(''.join(traceback.format_stack()))
                sys.exit(cd(f"assign_by_c12 allows either to consider\n\
                    * both\n\
                    * inner\n\
                    * outer\n\
                    '{location}' is not a valid choice."))


            if len(self.inner.get_data["id"])>=len(self.outer.get_data["id"]):
                randomizer=np.random.permutation(self.inner.get_data["id"])
            else:
                randomizer=np.random.permutation(self.outer.get_data["id"])

            for loc in locations:
                if loc.h_Name=="inner":
                    turn=-1
                else:
                    turn=1
                lipids={}
                for item in lines:
                    lipids[item[0]]=[round(float(item[2])*N,0)+1,float(item[3])]
                for index in randomizer:
                    domain=0
                    lipid_probabilities={}
                    Cs=np.asarray([loc.get_data["C1"][index],loc.get_data["C2"][index]])
                    for key in lipids:
                        Cs_input=np.asarray([lipids[key][1]])[0]
                        lipid_probabilities[key]=method

                    normalizer=np.sum(np.asarray(list(lipid_probabilities.values())))
                    if normalizer==0:
                        for key in lipid_probabilities:
                            lipid_probabilities[key]=1/len(lipid_probabilities.keys())
                            normalizer=1


                    for key in lipid_probabilities:
                        lipid_probabilities[key]=lipid_probabilities[key]/normalizer
                    domain=np.random.choice(list(lipid_probabilities.keys()),1,p=list(lipid_probabilities.values()))[0]
                    deleter=None
                    for key in lipids:
                        if int(domain)==int(key) and not unspecified_Number:
                            lipids[key][0]=lipids[key][0]-1
                            if lipids[key][0]==0 and len(lipids.keys())>1:
                                deleter=key
                    if deleter is not None:
                        del lipids[deleter]
                        deleter=None
                    loc.get_data["domain_id"][index]=domain

    def _backup_path(self,path):
        """
        _backup_path is a helper function to create sensible backups for the methods that actually write out files, or folders. It should
        not be called by itself. It checks if the path already exists and creates an alterantive path to save, which would not overwrite
        anything.

        :param path: a path (to which something might be saved)

        :return: returns a better path for backup creation
        """
        if path[-1]=="/":
            path=path[:-1]
        path_split=path.split("/")
        if os.path.exists(path):
            path_split[-1]="#"+path_split[-1]
            path="/".join(path_split)
            path=self._backup_path(path)
            return path
        else:
            return path

    def _cat_to_one(self,Object):
        """
        _cat_to_one is a helper function, which recombines the data in inner or outer to be writeable again to the original form. It should
        not be called separately.

        :param Object: Object is a class object for the inner or outer supclasses of the point class.

        :return: Returns the concatenated information ready for writing.
        """
        result=np.zeros((len(Object.get_data),len(Object.get_data["id"])))
        index=0
        for key in Object.get_data:
            result[index]=Object.get_data[key]
            index+=1
        return result.T

    def write_folder(self):
        """
        write_folder takes the point class information and writes it back to the point folder it read it from. The method overwrites the folder only
        after creating a backup of the old one.
        """
        backup=self._backup_path(self.path)
        if backup!=self.path:
            shutil.copytree(self.path,backup)
            print(f"Found Folder {self.path}. Backup created as {backup}")
        if self.path[-1]!="/":
            self.path=self.path+"/"
        savers={"InnerBM.dat":self.inner,"OuterBM.dat":self.outer,"ExcData.dat":self.exclusions,"IncData.dat":self.inclusions}
        for key in savers:
            if key in ["InnerBM.dat","OuterBM.dat"]:
                all_in_one=self._cat_to_one(savers[key])
                fmt="".join(['%10d','%5d']+['%10.3f']*4+['%8.3f']*11)
                header=f"< Point NoPoints       {all_in_one.shape[0]}>\n< id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2  >\n< {key[:5]} >"
                if "Outer" in key:
                    header=self.box+header
                np.savetxt(self.path+key,np.round(all_in_one,3),header=header,comments='',encoding="UTF8",fmt=fmt)
            elif key =="ExcData.dat":
                pass
            elif key =="IncData.dat":
                try:
                    all_in_one=self._cat_to_one(savers[key])
                    header=f"< Inclusion NoInc       {all_in_one.shape[0]}   >\n< id typeid pointid lx ly lz  >"
                    fmt="".join(['%12d']*3+['%8.3f']*3)
                    np.savetxt(self.path+key,np.round(all_in_one,3),header=header,comments='',encoding="UTF8",fmt=fmt)
                except KeyError:
                    pass

    def write_input_str(self,output_file="input_from_PUC.str",input_file=None,ts2cg_input_str=None):

        """
        write_input_str writes a input.str file that is readable by TS2CG to continue the workflow.

        :param output_file: (default="input.str") The file to which the new input.str should be written. A backup system
        is in place, which only overwrites an existing input.str after creating a backup
        :param input_file: (default=None) An input file of some form is required which matches the layout described in
        the point class description or for assign_by_c12 or assign_by_c0. If no input file is given, the input file
        assigned in the point class is used. However, if neither provided nor found in the class, an error will be raised.
        :param ts2cg_input_str: (default=None) The ts2cg_input_str has to have the form of an input.str used for TS2CG. If
        additional parts like inclusions or exclusions are given in the file, they will be reintroduced to the new input.str
        with the updated domains. If only the [Lipids List] section is used in the old input.str, the inclusion here is redundant.
        """
        if os.path.isdir(output_file):
            if output_file[:-1]=="/":
                output_file=output_file+"input_from_PUC.str"
            else:
                output_file=output_file+"/input_from_PUC.str"
        path=self._backup_path(output_file)
        if path!=output_file:
            shutil.copy(output_file,path)
            print(f"Found File {output_file}. Backup created as {path}")
        if input_file is not None:
            self.input_file=input_file
        elif self.input_file is None:
            print(traceback.format_exc())
            sys.exit(cd("The path to the input file has not been provided and not been set in the class, please supply\n\
            * the path with the input_file parameter."))
        with open(self.input_file,"r",encoding="UTF8") as f:
            lines=f.readlines()
            lines=[line for line in lines if line.strip()]
            lines=[line[:line.find(";")].split() for line in lines if line[0] != ";"]

        with open(output_file,"w",encoding="UTF8") as f:
            f.write("[Lipids List]\n")
            for line in lines:
                f.write(f"Domain {line[0]}\n")
                f.write(f"{line[1]} 1 1 {line[4]}\n")
                f.write("End\n")

        if ts2cg_input_str is not None:
            with open(ts2cg_input_str,"r",encoding="UTF8") as g:
                lines=g.readlines()

            with open(output_file,"a",encoding="UTF8") as f:
                block=False
                for line in lines:
                    if "Domain" in line:
                        block=True
                    elif "End" in line:
                        block=False
                    elif not block and "Lipids List" not in line:
                        f.write(line)

    def _wiggle(self,xyz,wiggle):
        random_perturbation=np.random.uniform(xyz-wiggle,xyz+wiggle,xyz.shape)
        xyz=xyz[0]+random_perturbation[0]

        return xyz/np.linalg.norm(xyz)


    def _corrected(self,to_check,reference):
        while len(to_check)<len(reference):
            to_check=to_check+[0]
        return to_check

    def set_protein_inclusion_number(self,number,Type=None,collision_Distance=0,domain=None, n_Vector=[0,0,0],wiggle=0):
        """
        Use this method to set the number of protein copies to a specified value in the Inclusions section.
        The id corresponds to the id in the outer membrane.

        :param number: The number of proteins there should be in the inclusions. Can also be a list.
        :param Type: (default=None) If a type is set, the number will correspond to the number of protein type as specified
        in the inclusion file. If no type is set, the number is split over all types according to the ratio of their
        appearance in the inclusion file.
        :param collision_Distance: (default=0) The placement will be restricted so a protein is not placed within
        this value to another protein. Additionally a list can be given to collision_Distance, which would contain
        the distance per protein type, i.e. [5,3,2] would say, that type 1 needs minimum space of 5 nm, 2 of 3, etc.
        :param domain: (default=None) If a domain is set, the inclusions will only be placed in that domain. Here,
        also a list of domains is possible to allow placement in multiple domains.
        :param n_Vector: (default=[0,0,0]) The relative position of the protein type is set. Set n_Vector to None
        to take the n_Vector from the inclusion file.
        :param wiggle: (default=0) The relative position is wiggled a little per inclusion.
        """

        #Make everything a list to unify the code disregarding if the paramters are lists or not
        if not isinstance(Type,list):
            Type=[Type]
        if not isinstance(number,list):
            number=[number]
        if not isinstance(collision_Distance,list):
            collision_Distance=[collision_Distance]
        if not isinstance(domain,list):
            domain=[domain]
        if n_Vector is not None:
            if not any(isinstance(item,list) for item in n_Vector):
                n_Vector=[n_Vector]
        if not isinstance(wiggle,list):
            wiggle=[wiggle]
        number=self._corrected(number,Type)
        collision_Distance=self._corrected(collision_Distance,Type)
        wiggle=self._corrected(wiggle,Type)

        #baseline items will have the form [typeid,desired number,averaged position,the allowed wiggle]
        baseline=np.zeros((len(Type),7))
        ratio=[]
        inc_baseline=[]
        for i,item in enumerate(Type):
            baseline[i][0]=item
            ratio.append(self.inclusions.raw.T[self.inclusions.raw.T[:,1]==item])
            mean=np.mean(self.inclusions.raw.T[self.inclusions.raw.T[:,1]==item][:,3:],axis=0)
            baseline[i][4:]=mean
            baseline[i][1]=number[i]
            baseline[i][2]=wiggle[i]
            baseline[i][3]=collision_Distance[i]

            inc_baseline=inc_baseline+([item]*number[i])


        outer_baseline=self.outer.raw.T[:,[0,1,3,4,5]]
        outer_baseline=outer_baseline[np.isin(outer_baseline[:,1],domain)]
        outer_baseline=outer_baseline[np.random.permutation(outer_baseline.shape[0])]
        inc_baseline=np.asarray(inc_baseline)[np.random.permutation(len(inc_baseline))]


        collisions=[]
        colliding_distances=[]
        new_raw=[]
        index=0

        for item in outer_baseline:
            placeable=True
            for i in range(len(collisions)):
                if np.linalg.norm(item[2:]-collisions[i])<(baseline[baseline[:,0]==inc_baseline[index]][:,3]+colliding_distances[i])/2:
                    placeable=False
                    break
            if placeable:
                collisions.append(item[2:])
                colliding_distances.append(baseline[baseline[:,0]==inc_baseline[index]][:,3])
                wiggled=self._wiggle(xyz=baseline[baseline[:,0]==inc_baseline[index]][:,4:],wiggle=baseline[baseline[:,0]==inc_baseline[index]][:,2])
                new_raw.append([index,inc_baseline[index],item[0],wiggled[0],wiggled[1],wiggled[2]])
                index+=1
            if index==len(inc_baseline):
                break

        self.inclusions.make_data(np.asarray(new_raw).T)

    def _protein_from_input(self,Config):
        future_kwargs={"number":None,"Type":None,"collision_Distance":0,"domain":None,"n_Vector":[0,0,0],"wiggle":0}
        with open(Config,"r",encoding="UTF8") as f:
            lines=f.readlines()
            lines=[line for line in lines if line.strip()]
            lines=[line[:line.find(";")].split() for line in lines if line[0] != ";"]

            for key in future_kwargs:
                for item in lines:
                    if key in item:
                        future_kwargs[key]=line[line.find(":")+1:].strip().split(",")
                        for i,value in enumerate(future_kwargs[key]):
                            try:
                                future_kwargs[key][i]=int(value)
                            except ValueError:
                                future_kwargs[key][i]=float(value)
        if future_kwargs["number"] is None:
            print(''.join(traceback.format_stack()))
            sys.exit(cd(f"The minimum requirement is supplying an amount for a protein type"))

    def _find_dummies(self,dummys):
        dummys=dummys.split(",")
        for i,index in enumerate(dummys):
            try:
                dummys[i]=int(index)
            except ValueError:
                print("WARNING: a dummy protein could not be parsed.")
        print(dummys)
        return np.asarray(dummys)


    def domain_around_inclusion(self,radius,Type,domain=1,layer="both",dummy_Prot=""):
        Types=self.inclusions.raw[1]
        pointids=self.inclusions.raw[2]
        interest=np.asarray(pointids[Types==int(Type)],dtype=np.int32)
        print(interest)
        interest=np.concatenate((interest,self._find_dummies(dummy_Prot)))
        print(interest)
        interest_coords=[]
        for Id in interest:
            interest_coords.append(self.outer.raw.T[int(Id),3:6])

        if layer=="both":
            raw_trans=self.outer.raw.T
            raw_trans_inner=self.inner.raw.T
        elif layer=="inner":
            raw_trans=self.inner.raw.T
        else:
            raw_trans=self.outer.raw.T

        for location in interest_coords:
            for i,line in enumerate(raw_trans):
                coords=line[3:6]
                if np.linalg.norm(coords-location)<float(radius):
                    print("domain changed")
                    raw_trans[i,1]=domain
                    if layer=="both":
                        raw_trans_inner[i,1]=domain

        if layer=="both":
            self.outer.raw=raw_trans.T
            self.inner.raw=raw_trans_inner.T
            self.inner.make_dict()
            self.outer.make_dict()
        elif layer=="inner":
            self.inner.raw=raw_trans.T
            self.inner.make_dict()
        else:
            self.outer.raw=raw_trans.T
            self.outer.make_dict()


    @classmethod
    def DOP(cls, args):
        parser=argparse.ArgumentParser()
        parser = argparse.ArgumentParser(description="A tool to directly alter the domains to place lipids according to a preferred curvature",formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument('-i','--input',type=str,default="domain.txt",help="""A path to the domain.txt input file. A file could for instance, look like this:
            ; domain lipid percentage c0 density
            0 POPC .3 0.179 0.64
            1 POPD .4 0.613 0.64
            2 POPG .3 0.629 0.64""")
        parser.add_argument('-k','--k',default=1,help="Sets the factor in the exponent for the c0 approach")
        parser.add_argument('-p','--path',default="point/",help="Specify the path to the point folder")
        parser.add_argument('-l','--location',default='both',help="Choose which monolayer is altered: both, upper, lower")
        parser.add_argument('-ni','--new_TS2CG',default=None,help="Path to write a new TS2CG input file.")
        parser.add_argument('-oi','--old_TS2CG',default=None,help="Supply a path to the TS2CG input.str")
        parser.add_argument('-I','--Ignore_lipid_number',action='store_true',default=False,help="Set to ignore the number of lipids given in the input file.")
        parser.add_argument('-n','--new',action='store_true',help='Set new to initialize a new point folder from scratch, will fail if the folder already exists')

        args=parser.parse_args(args)

        PointFolder=cls(path=args.path,k=args.k,new=args.new)
        PointFolder.assign_by_c0(args.input,location=args.location,unspecified_Number=args.Ignore_lipid_number)
        if args.new_TS2CG is not None:
            PointFolder.write_input_str(output_file=args.new_TS2CG,input_file=args.input,ts2cg_input_str=args.old_TS2CG)
        PointFolder.write_folder()

    @classmethod
    def INU(cls, args):
        parser=argparse.ArgumentParser()
        parser = argparse.ArgumentParser(description="A tool to directly alter the number of protein inclusions in the membrane",formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument('-i','--input',type=str,default="proteins.txt",help="""A path to the proteins.txt input file. The file needs to contain
        the parameters for the protein placement. The input file should read like:
            ; this is a comment
            number: 1,2,3
            Type: 3,4  ; default is None
            collision_Distance: 4.5,6.2 ;default is 0
            domain: 4,5 ; default is None
            n_Vector: 1,3,4 ; xyz coordinates of a vector that will be normalized later. Default is 0,0,0
            wiggle: 11.2 ; A value that wiggles in a normal distribution around n_vector for some variance in protein orientation. Default is 0
            """)
        parser.add_argument('-p','--path',default="point/",help="Specify the path to the point folder")
        parser.add_argument('-ni','--new_TS2CG',default=None,help="Path to write a new TS2CG input file.")
        parser.add_argument('-oi','--old_TS2CG',default=None,help="Supply a path to the TS2CG input.str")

        args=parser.parse_args(args)
       
        PointFolder=cls(path=args.path)
        input_args=PointFolder.protein_from_input(args.input)
        PointFolder._set_protein_inclusion_number(self,**input_args)
        PointFolder.protein_altered=True
        if args.new_TS2CG is not None:
            PointFolder.write_input_str(output_file=args.new_TS2CG,input_file=args.input,ts2cg_input_str=args.old_TS2CG)
        PointFolder.write_folder()

    @classmethod
    def DAI(cls, args):
        parser=argparse.ArgumentParser()
        parser = argparse.ArgumentParser(description="A tool to directly alter the domain composition circular around inclusions.",formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument('-p','--path',default="point/",help="Specify the path to the point folder")
        parser.add_argument('-ni','--new_TS2CG',default=None,help="Path to write a new TS2CG input file.")
        parser.add_argument('-oi','--old_TS2CG',default=None,help="Supply a path to the TS2CG input.str")
        parser.add_argument('-r','--radius',default=1,type=float,help="The radius around a protein in which domains should be changed.")
        parser.add_argument('-T','--Type',default=1,type=int,help="The protein type around which domains should be changed.")
        parser.add_argument('-d','--Domain',default=1,type=int,help="The domain number that should be set around the protein.")
        parser.add_argument('-L','--Layer',default="both",help="Choose which membrane layer to alter. Default is both")
        parser.add_argument('-dummy','--dummy',default="",help="Crate a dummy protein to place a circular domain around it. Excepts pointids like 3,7,22")


        args=parser.parse_args(args)
       
        PointFolder=cls(path=args.path)
        input_args=PointFolder.domain_around_inclusion(args.radius,args.Type,args.Domain,args.Layer,args.dummy)
        if args.new_TS2CG is not None:
            PointFolder.write_input_str(output_file=args.new_TS2CG,input_file=args.input,ts2cg_input_str=args.old_TS2CG)
        PointFolder.write_folder()
        

