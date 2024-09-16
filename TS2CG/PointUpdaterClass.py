import TS2CG.PointUpdater.point as p
import argparse
from inspect import cleandoc as cd


def PUC(args):
    """
    PUC is short for Point Updater Class and is the most straight forward way to use the PointUpdater.
    PUC is invoked like the other steps in TS2CG by 'TS2CG PUC --flags', where the flags specify the
    needed input.
    """
    parser=argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description="A tool to directly alter the domains to place lipids according to a preferred curvature",formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i','--input',type=str,default="domain.txt",help="""A path to the domain_input.txt input file. A file could for instance, look like this:
        ; domain lipid percentage c0 density
        0 POPC .3 0.179 0.64
        1 POPD .4 0.613 0.64
        2 POPG .3 0.629 0.64""")
    parser.add_argument('-c','--c12',action='store_true',help="Use the c12 approach to set the domains instead of the default c0 approach.")
    parser.add_argument('-p','--path',default="point/",help="Specify the path to the point folder")
    parser.add_argument('-l','--location',default='both',help="Choose which monolayer is altered: both, upper, lower")
    parser.add_argument('-ni','--new_TS2CG',default=None,help="Path to write a new TS2CG input file.")
    parser.add_argument('-oi','--old_TS2CG',default=None,help="Supply a path to the TS2CG input.str")

    args=parser.parse_args(args)
    
    PointFolder=p.point(args.path)
    if args.c12:
        PointFolder.assign_by_c12(args.input,location=args.location)
    else:
        PointFolder.assign_by_c0(args.input,location=args.location)
    if args.new_TS2CG is not None:
        PointFolder.write_input_str(output_file=args.new_TS2CG,input_file=args.input,ts2cg_input_str=args.old_TS2CG)
    PointFolder.write_folder()
    
    





if __name__=="__main__":
    pass
    #PUC("-i PointUpdater/domain_input.txt -p ../../point_python/point/")