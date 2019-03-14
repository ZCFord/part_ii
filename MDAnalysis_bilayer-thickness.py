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

parser = argparse.ArgumentParser('Calculates the bilayer thickness and average z boundaries')

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
	downP = L.groups(0)
	upP = L.groups(1)
	upPxyz = upP.positions
	upPz = upPxyz[:,2] #selects z coordinates
	upPz_av = np.mean(upPz)
	downPxyz = downP.positions
	downPz = downPxyz[:,2] 
	downPz_av = np.mean(downPz)
	bilayer_thickness = downPz_av - upPz_av #difference in average upper and lower z coordinates
	return (bilayer_thickness)

def z_upper (A) : #finds the upper z boundary of the membrane
	upP = L.groups(1)
	upPxyz = upP.positions
	upPz = upPxyz[:,2]
	upPz_av = np.mean(upPz)
	return (upPz_av)

def z_lower (A) : #finds the upper z boundary of the membrane
	downP = L.groups(0)
	downPxyz = downP.positions
	downPz = downPxyz[:,2]
	downPz_av = np.mean(downPz)
	return (downPz_av)

def biup_std (A) : #standard deviation of upper leaflet
	upP = L.groups(1)
	upPxyz = upP.positions
	upPz = upPxyz[:,2]
	upPz_std = np.std(upPz)
	return (upPz_std)

def bidown_std (A) : #standard deviation of lower leaflet
	downP = L.groups(0)
	downPxyz = downP.positions
	downPz = downPxyz[:,2]
	downPz_std = np.std(downPz)
	return (downPz_std)

def total_std (A, B) : #combined standard deviations
	x2 = A**2 + B**2
	x = math.sqrt(x2)
	return x

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

##########
# CALCULATIONS
##########

bilayer = []
high_std = []
low_std = []

u = MDAnalysis.Universe(args.f, args.t)
L = MDAnalysis.analysis.leaflet.LeafletFinder(u, 'name P*') #selects P atoms and puts in two leaflet groups

for ts in u.trajectory :
	bilayer.append(bilayer_thickness(u))
	up = biup_std(u)
	down = bidown_std(u)
	tot = total_std(up, down) #up and down are A and B
	high = bilayer_thickness(u) + tot
	low = bilayer_thickness(u) - tot
	high_std.append(high)
	low_std.append(low)

av_run = running_mean(bilayer, 10)
high = running_mean(high_std, 10)
low = running_mean(low_std, 10)

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
x = np.arange(4,96,1)
plt.plot(bilayer, label = 'Bilayer Thickness', color = 'blue', linewidth = 3)
plt.plot(x, av_run, label = '10ns Running Average', color = 'blue', linewidth = 2, linestyle = 'dashed')
plt.plot(x, high, label = '+1 Standard Deviation', color = 'green', linewidth = 1)
plt.plot(x, low, label = '-1 Standard Deviation', color = 'red', linewidth = 1)
v = [0,100,20,50]
plt.axis(v)
plt.fill_between(x, high, low, facecolor = 'blue', alpha = 0.2)
plt.axhline(y=(bilayer[0]), color = 'black', linestyle = 'dashed', linewidth = 1, label = 'Starting Thickness')
plt.ylabel("Bilayer Thickness ($\AA$)", fontsize = '24')
plt.xlabel("Time (ns)", fontsize = '24')
plt.legend(fontsize = 'x-large', loc = 'best')
plt.savefig("{}_bilayer.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_bilayer.svg".format(args.o), format='svg', dpi=300)                               
