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
	help = "First structure (pdb)"
)

parser.add_argument(
	"-ft",
	type = str,	
	help = "First trajectory (xtc)"
)

parser.add_argument(
	"-fl",
	type = str,
	help = "First label"
)

parser.add_argument(
	"-s",
	type = str,
	help = "Second structure (pdb)"
)

parser.add_argument(
	"-st",
	type = str,	
	help = "Second trajectory (xtc)"
)

parser.add_argument(
	"-sl",
	type = str,
	help = "Second label"
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

def bilayer_thickness (A, B) : # finds the bilayer thickness from a file or universe input
	downP = B.groups(0)
	upP = B.groups(1)
	upPxyz = upP.coordinates()
	upPz = upPxyz[:,2] #selects z coordinates
	upPz_av = np.mean(upPz)
	downPxyz = downP.coordinates()
	downPz = downPxyz[:,2] 
	downPz_av = np.mean(downPz)
	bilayer_thickness = downPz_av - upPz_av #difference in average upper and lower z coordinates
	return (bilayer_thickness)

def z_upper (A, B) : #finds the upper z boundary of the membrane
	upP = B.groups(1)
	upPxyz = upP.coordinates()
	upPz = upPxyz[:,2]
	upPz_av = np.mean(upPz)
	return (upPz_av)

def z_lower (A, B) : #finds the upper z boundary of the membrane
	downP = B.groups(0)
	downPxyz = downP.coordinates()
	downPz = downPxyz[:,2]
	downPz_av = np.mean(downPz)
	return (downPz_av)

def biup_std (A, B) : #standard deviation of upper leaflet
	upP = B.groups(1)
	upPxyz = upP.coordinates()
	upPz = upPxyz[:,2]
	upPz_std = np.std(upPz)
	return (upPz_std)

def bidown_std (A, B) : #standard deviation of lower leaflet
	downP = B.groups(0)
	downPxyz = downP.coordinates()
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

firstbilayer = []
high_std_1 = []
low_std_1 = []

u = MDAnalysis.Universe(args.f, args.ft)
L = MDAnalysis.analysis.leaflet.LeafletFinder(u, 'name P*') #selects P atoms and puts in two leaflet groups

for ts in u.trajectory :
	firstbilayer.append(bilayer_thickness(u, L))
	up = biup_std(u, L)
	down = bidown_std(u, L)
	tot = total_std(up, down) #up and down are A and B
	high = bilayer_thickness(u, L) + tot
	low = bilayer_thickness(u, L) - tot
	high_std_1.append(high)
	low_std_1.append(low)

first_av_run = running_mean(firstbilayer, 10)
high_1 = running_mean(high_std_1, 10)
low_1 = running_mean(low_std_1, 10)

secondbilayer = []
high_std_2 = []
low_std_2 = []

v = MDAnalysis.Universe(args.s, args.st)
J = MDAnalysis.analysis.leaflet.LeafletFinder(v, 'name P*') #selects P atoms and puts in two leaflet groups

for ts in v.trajectory :
	secondbilayer.append(bilayer_thickness(v, J))
	up = biup_std(v, J)
	down = bidown_std(v, J)
	tot = total_std(up, down) #up and down are A and B
	high = bilayer_thickness(v, J) + tot
	low = bilayer_thickness(v, J) - tot
	high_std_2.append(high)
	low_std_2.append(low)

second_av_run = running_mean(secondbilayer, 10)
high_2 = running_mean(high_std_2, 10)
low_2 = running_mean(low_std_2, 10)

##########
# MEAN BILAYER BOUNDARIES
##########

upper_1 = []
lower_1 = []

for ts in u.trajectory :
	upper_1.append(z_upper(u, L))
	lower_1.append(z_lower(u, L))

up_av_1 = np.mean(upper_1)
lo_av_1 = np.mean(lower_1)
print(up_av_1)
print(lo_av_1)

upper_2 = []
lower_2 = []

for ts in v.trajectory :
	upper_2.append(z_upper(v, J))
	lower_2.append(z_lower(v, J))

up_av_2 = np.mean(upper_2)
lo_av_2 = np.mean(lower_2)
print(up_av_2)
print(lo_av_2)

##########
# PLOTTING
##########

plt.figure(figsize=(16,8))
x = np.arange(4,96,1)
plt.plot(firstbilayer, label = '{} Bilayer Thickness'.format(args.fl), color = 'blue', linewidth = 3)
plt.plot(x, first_av_run, label = '10ns Running Average', color = 'blue', linewidth = 2, linestyle = 'dashed')
plt.plot(secondbilayer, label = '{} Bilayer Thickness'.format(args.sl), color = 'red', linewidth = 3)
plt.plot(x, second_av_run, label = '10ns Running Average', color = 'red', linewidth = 2, linestyle = 'dashed')
#plt.plot(x, high, label = '+1 Standard Deviation', color = 'green', linewidth = 1)
#plt.plot(x, low, label = '-1 Standard Deviation', color = 'red', linewidth = 1)
v = [0,100,20,50]
plt.axis(v)
#plt.fill_between(x, high, low, facecolor = 'blue', alpha = 0.2)
#plt.axhline(y=(bilayer[0]), color = 'black', linestyle = 'dashed', linewidth = 1, label = 'Starting Thickness')
plt.ylabel("Bilayer Thickness ($\AA$)", fontsize = '24')
plt.xlabel("Time (ns)", fontsize = '24')
plt.legend(fontsize = 'x-large', loc = 'best')
plt.savefig("{}_bilayer.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_bilayer.svg".format(args.o), format='svg', dpi=300)                               
