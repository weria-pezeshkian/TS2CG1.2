import numpy as np
from TS2CG.PointUpdater import point as p

class Init_Point(p.point):
    def __init__(self):
        pass


def builder(location="both"):
    p=Init_Point()
    if location in ["inner","outer"]:
        pass
    else:
        inner=p.InnerOrOuter(init=True,h_Name="inner")