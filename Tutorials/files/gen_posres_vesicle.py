import sys
import numpy as np
import MDAnalysis as mda

input = sys.argv[1]
output = sys.argv[2]

u = mda.Universe(input)
cog = u.atoms.center_of_geometry()
u.atoms.positions = cog
u.atoms.write(output, reindex=False)
