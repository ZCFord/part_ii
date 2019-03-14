import MDAnalysis
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import argparse

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False

##########
# ARGUMENT HANDLING
##########

parser = argparse.ArgumentParser('Calculates the z boundaries of the selectivity filter')

parser.add_argument(
	"-f",
	type = str,
	default = './mdord.pdb',
	help = "Structure (pdb)"
)

parser.add_argument(
	"-t",
	type = str,	
	default = 'md-c.xtc',
	help = "Trajectory (xtc)"
)

parser.add_argument(
	"-o",
	type = str,
	default = None,
	help = "output name"
)

args = parser.parse_args()

##########
# CALCULATIONS
##########

## want to get average z position of Ts and upper Gs in selectivity filter

u = MDAnalysis.Universe(args.f, args.t)
protein = u.select_atoms('protein')

gly = []
thr = []

for ts in u.trajectory :
	glys = protein.select_atoms('resname GLY and resid 176 or resid 285').coordinates()
	gly.append(np.mean(glys[:,2]))
	thrs = protein.select_atoms('resname THR and resid 172 or resid 281').coordinates()
	thr.append(np.mean(thrs[:,2]))

print np.mean(gly)
print np.mean(thr)

plt.plot(gly)
plt.plot(thr)

