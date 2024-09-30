import numpy as np
import sys
import traceback
from inspect import cleandoc as cd
import os,shutil


class point():
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
        def __init__(self,data):
            """
            :param data: The data of the appropriate dat file as has been read by read_BM()
            """
            names="id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2".split()
            self.get_data={}
            self.raw=data
            for i,row in enumerate(self.raw):
                self.get_data[names[i]]=row
            for name in ["id","domain_id"]:
                self.get_data[name]=np.asarray(self.get_data[name],dtype=int)
    
    class Inclusions():
        """
        Represents the IncData.dat from the point folder.
        """
        def __init__(self,data):
            """
            :param data: The data of the appropriate dat file as has been read by read_BM()
            """
            names="id typeid pointid lx ly lz".split()
            self.get_data={}
            self.raw=data
            try:
                for i,row in enumerate(self.raw):
                    self.get_data[names[i]]=row
            except TypeError:
                pass

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


    def __init__(self, path="point/"):
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
        or
            ; domain lipid percentage c1 c2 density
            0 POPC .3 0.179 0.150 0.64
            1 POPD .4 0.613 0.481 0.64
            2 POPG .3 0.629 0.472 0.64
        depending on the alteration planned later.
        If necessarry k can be adjusted to correct the exponent in the formula, if c0 is used to sort lipids by curvature.

        :param path: (default="point/)") gives the path to the point folder to be read out.
        """
        self.box="Box XXX XXX XXX"
        if path[-1]=="/":
            self.path=path
        else:
            self.path=path+"/"
        self.input_file=None
        try:
            self.inner=self.InnerOrOuter(self.read_BM(self.path+"InnerBM.dat"))
            self.outer=self.InnerOrOuter(self.read_BM(self.path+"OuterBM.dat"))
            self.inclusions=self.Inclusions(self.read_BM(self.path+"IncData.dat",coords=False))
            self.exclusions=self.Exclusions(self.read_BM(self.path+"ExcData.dat",coords=False))
        except:
            print(traceback.format_exc())
            sys.exit(cd("Something went wrong initializing the point class. Please make sure that\n\
            * path is set to the point folder\n\
            * the point folder contains the files InnerBM.dat, OuterBM.dat, IncData.dat, and ExcData.dat"))
        self.k=1

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
                lipids={}
                for item in lines:
                    lipids[item[0]]=[round(float(item[2])*N,0)+1,float(item[3])]
                for index in randomizer:
                    domain=0
                    lipid_probabilities={}
                    Cs=np.asarray([loc.get_data["C1"][index],loc.get_data["C2"][index]])
                    for key in lipids:
                        Cs_input=np.asarray([lipids[key][1]])[0]
                        lipid_probabilities[key]=np.exp(-self.k*(Cs[0]+Cs[1]-Cs_input)**2)
                    
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

    def assign_by_c12(self,input_file,location="both",unspecified_Number=False):
        """
        the function assign_by_c12 reads an input file of form (axample)
            ; domain lipid percentage c1 c2 density
            0 POPC .3 0.179 0.150 0.64
            1 POPD .4 0.613 0.481 0.64
            2 POPG .3 0.629 0.472 0.64
        which would set domains 0, 1, 2 for the lipids randomly with their corresponding percentage such that the given c1, c2 values are
        closest to the c1, c2 values of the membrane in an Euclidean distance for the yet to be occupied spaces.

        :param input_file: path to the input file
        :param location: (default="both") allows the alteration of just the *upper*,*lower*, or *both* monolayer.
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
                lipids={}
                for item in lines:
                    lipids[item[0]]=[round(float(item[2])*N,0),float(item[3]),float(item[4])]
                for index in randomizer:
                    domain=0
                    distance=np.inf
                    Cs=np.asarray([loc.get_data["C1"][index],loc.get_data["C2"][index]])
                    for key in lipids:
                        Cs_input=np.asarray([lipids[key][1:]])
                        if lipids[key][0]!=0 and np.linalg.norm(Cs-Cs_input)<distance: #find a smarter idea than cartesian distance for curvature
                            distance=np.linalg.norm(Cs-Cs_input)
                            domain=int(key)
                    for key in lipids:
                        if domain==int(key) and not unspecified_Number:
                            lipids[key][0]=lipids[key][0]-1
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
                np.savetxt(self.path+key,np.round(all_in_one,3),header=header,encoding="UTF8",fmt=fmt)
            elif key =="ExcData.dat":
                pass
                #TODO: need to find out what to do with the ExcData
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


if __name__=="__main__":
    pass
    #a=point(path="../../../point_python/point/")
    #a.write_input_str(output_file="input.str",input_file="domain_input.txt",ts2cg_input_str="../../../MobiusMembrane/input.str")
    #a.assign_domains_by_C("domain_input.txt",location="both",c0=True)
    #a.write_folder()

