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

parser = argparse.ArgumentParser('Plotting LCA')

parser.add_argument(
	"-n",
	type = str,
	help = "not stretched"
)

parser.add_argument(
	"-s",
	type = str,	
	help = "stretched"
)

parser.add_argument(
	"-o",
	type = str,	
	help = "output name (LIPID)"
)

args = parser.parse_args()

stretch = np.loadtxt(args.s)
nostretch = np.loadtxt(args.n)

nstretch = stretch[:,0]-(stretch[-1,0]/2)
NOstretch = nostretch[:,0]-(nostretch[-1,0]/2)

plt.figure(figsize=(10,8))
plt.title('Bilayer Pressure Profile', fontsize = 30)
plt.plot(NOstretch*10, (nostretch[:,1]+nostretch[:,5])/2, label = 'No Stretch', color = 'blue')
plt.plot(nstretch*10, (stretch[:,1]+stretch[:,5])/2, label = '-40 bar Stretch', color = 'red')
v = [-50,50,-500,800]
plt.axis(v)
plt.ylabel("P (bar)", fontsize = '24')
plt.xlabel("Z ($\AA$)", fontsize = '24')
plt.legend(fontsize = 'x-large', loc = 'best')
plt.savefig("{}stress_comp.png".format(args.o), format='png', dpi=300)
plt.savefig("{}stress_comp.svg".format(args.o), format='svg', dpi=300)   

