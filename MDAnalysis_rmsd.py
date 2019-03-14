#### NOTE: THIS IS NOT FOR THE CHIMERA

import MDAnalysis
import MDAnalysis.analysis.rms
import matplotlib
import matplotlib.pyplot as plt
import argparse
from matplotlib import rc, rcParams

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False

##########
# ARGUMENT HANDLING
##########

parser = argparse.ArgumentParser('Calculates RMSD between two structures')

parser.add_argument(
	"-f",
	type = str,
	required=True,
	help="first input; if only using one trajectory it should be for this input"
)

parser.add_argument(
	"-s",
	type = str,	
	required=True,
	help="second input"
)

parser.add_argument(
	"-ft",
	type = str,
	help="first trajectory"
)

parser.add_argument(
	"-o",
	type = str,
	default = None,
	help="output name"
)

#parser.add_argument(
#	"-st",
#	type = str,
#	help="second trajectory; requires ft input as well"
#)

args = parser.parse_args()

#########
# MAKE UNIVERSE
#########

print ('Have you remembered to fix the rotation of your protein? Remember to use mdord!')

# first structure
if args.ft == None :
	first = MDAnalysis.Universe(args.f)
else :
	first = MDAnalysis.Universe(args.f, args.ft)

# second structure
#if args.st == None :
#	second = MDAnalysis.Universe(args.s)
#else :
#	second = MDAnalysis.Universe(args.s, args.st)

second = MDAnalysis.Universe(args.s)

#########
# ONE FRAME
#########

if args.ft == None : # and args.st == None :
	protein_1 = first.select_atoms('protein')
	protein_2 = second.select_atoms('protein')
	sel_1 = protein_1.select_atoms('(resid 72-98 or resid 183-229 or resid 232-264 or resid 299-330) and name CA')
	sel_2 = protein_2.select_atoms('(resid 72-98 or resid 183-229 or resid 232-264 or resid 299-330) and name CA')
	one = sel_1.positions.copy()
	two = sel_2.positions.copy()
	print(MDAnalysis.analysis.rms.rmsd(one, two, superposition=True))
	print("ay up")

#########
# ONE TRAJECTORY, ONE STATIC STRUCTURE
#########

else :
	protein_2 = second.select_atoms('protein')
	sel_2 = protein_2.select_atoms('(resid 72-98 or resid 183-229 or resid 232-264 or resid 299-330) and name CA')
	two = sel_2.positions.copy()
	print("wotcher")
	
	RMSD_traj = []
	for ts in first.trajectory :
		protein_1 = first.select_atoms('protein')
		sel_1 = protein_1.select_atoms('(resid 72-98 or resid 183-229 or resid 232-264 or resid 299-330) and name CA')
		one = sel_1.positions.copy()
		RMSD_traj.append(MDAnalysis.analysis.rms.rmsd(one,two, superposition=True))

	plt.figure(figsize=(10,8))
	plt.plot(RMSD_traj, label = 'RMSD', color = 'blue', linewidth = 2)
	v = [0,100,1,5]
	plt.axis(v)
	plt.title("RMSD of a trajectory compared to a static structure", fontsize = 24)
	plt.ylabel("RMSD ($\AA$)", fontsize = '24')
	plt.xlabel("Time (ns)", fontsize = '24')
	plt.savefig("{}_rmsd.png".format(args.o), format='png', dpi=300)
	plt.savefig("{}_rmsd.svg".format(args.o), format='svg', dpi=300)

#########
# TWO TRAJECTORIES
#########

#else :
#	print("oy oy")
#	one = []
#	two = []	
#	for ts in first.trajectory :
#		protein_1 = first.select_atoms('protein')
#		sel_1 = protein_1.select_atoms('resid 1-250')
#		one.append(sel_1.positions.copy())
#
#	for ts in second.trajectory :
#		protein_2 = second.select_atoms('protein')
#		sel_2 = protein_2.select_atoms('resid 1-250')
#		two.append(sel_2.positions.copy())
#
#	RMSD = MDAnalysis.analysis.rms.rmsd(one, two)

#	plt.figure(figsize=(10,8))
#	plt.plot(RMSD, label = 'RMSD', color = 'blue', linewidth = 2)
#	v = [0,100,30,50]
#	plt.axis(v)
#	plt.title("RMSD of two trajectories", fontsize = 36)
#	plt.ylabel("RMSD ($\AA$)", fontsize = '24')
#	plt.xlabel("Time (ns)", fontsize = '24')
#	plt.savefig("rmsd.png", format='png', dpi=300)
#	plt.savefig("rmsd.svg", format='svg', dpi=300)                               
