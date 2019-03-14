import MDAnalysis
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis.analysis.leaflet
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

parser = argparse.ArgumentParser('Calculates angle of M2, M3, and M4 helices')

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
# DEFINITIONS
##########

def bilayer_thickness (A) : # finds the bilayer thickness from a file or universe input
	L = MDAnalysis.analysis.leaflet.LeafletFinder(A, 'name P*') #selects P atoms and puts in two leaflet groups
	downP = L.groups(0)
	upP = L.groups(1)
	upPxyz = upP.coordinates()
	upPz = upPxyz[:,2] #selects z coordinates
	upPz_av = np.mean(upPz)
	downPxyz = downP.coordinates()
	downPz = downPxyz[:,2] 
	downPz_av = np.mean(downPz)
	bilayer_thickness = downPz_av - upPz_av #difference in average upper and lower z coordinates
	return (bilayer_thickness)

def z_upper (A) : #finds the upper z boundary of the membrane
	L = MDAnalysis.analysis.leaflet.LeafletFinder(A, 'name P*')
	upP = L.groups(1)
	upPxyz = upP.coordinates()
	upPz = upPxyz[:,2]
	upPz_av = np.mean(upPz)
	return (upPz_av)

def z_lower (A) : #finds the upper z boundary of the membrane
	L = MDAnalysis.analysis.leaflet.LeafletFinder(A, 'name P*')
	downP = L.groups(0)
	downPxyz = downP.coordinates()
	downPz = downPxyz[:,2]
	downPz_av = np.mean(downPz)
	return (downPz_av)

def biup_std (A) : #standard deviation of upper leaflet
	L = MDAnalysis.analysis.leaflet.LeafletFinder(A, 'name P*')
	upP = L.groups(1)
	upPxyz = upP.coordinates()
	upPz = upPxyz[:,2]
	upPz_std = np.std(upPz)
	return (upPz_std)

def bidown_std (A) : #standard deviation of lower leaflet
	L = MDAnalysis.analysis.leaflet.LeafletFinder(A, 'name P*')
	downP = L.groups(0)
	downPxyz = downP.coordinates()
	downPz = downPxyz[:,2]
	downPz_std = np.std(downPz)
	return (downPz_std)

def total_std (A, B) : #combined standard deviations
	x2 = A**2 + B**2
	x = math.sqrt(x2)
	return x

##########
# CALCULATIONS
##########

bilayer = []
high_std = []
low_std = []

u = MDAnalysis.Universe(args.f, args.t)

for ts in u.trajectory :
	bilayer.append(bilayer_thickness(u))
	up = biup_std(u)
	down = bidown_std(u)
	tot = total_std(up, down) #up and down are A and B
	high = bilayer_thickness(u) + tot
	low = bilayer_thickness(u) - tot
	high_std.append(high)
	low_std.append(low)

##########
# MEAN BILAYER BOUNDARIES
##########

upper = []
lower = []

for ts in u.trajectory :
	upper.append(z_upper(u))
	lower.append(z_lower(u))

up_av = np.mean(upper)
lo_av = np.mean(lower)
print(up_av)
print(lo_av)

##########
# PLOTTING
##########

plt.figure(figsize=(10,8))
plt.plot(bilayer, label = 'Bilayer Thickness', color = 'blue', linewidth = 3)
plt.plot(high_std, label = '+1 std', color = 'green', linewidth = 1)
plt.plot(low_std, label = '-1 std', color = 'red', linewidth = 1)
v = [0,100,20,50]
plt.axis(v)
x = np.arange(0,101,1)
plt.fill_between(x, high_std, low_std, facecolor = 'blue', alpha = 0.2)
plt.ylabel("Bilayer Thickness ($\AA$)", fontsize = '24')
plt.xlabel("Time (ns)", fontsize = '24')
plt.savefig("{}_bilayer.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_bilayer.svg".format(args.o), format='svg', dpi=300)                               
