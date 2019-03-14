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

parser = argparse.ArgumentParser('Plot RDF')

parser.add_argument(
	"-o",
	type = str,
	default = None,
	help = "output name"
)

args = parser.parse_args()

h_RDF = np.loadtxt('./headrdf.xvg', skiprows = 25)
t_RDF = np.loadtxt('./tailrdf.xvg', skiprows = 25)

plt.figure(figsize=(10,8))
plt.title('Radial Distribution', fontsize = 36)
plt.plot(h_RDF[:,0]*10, h_RDF[:,1], label = 'Lipid Headgroups', color = 'red', linewidth = 1)
plt.plot(h_RDF[:,0]*10, t_RDF[:,1], label = 'Lipid Tails', color = 'blue', linewidth = 1)
v = [0,14,0,35]
plt.axis(v)
plt.ylabel("RDF g(r)", fontsize = '24')
plt.xlabel("Radius ($\AA$)", fontsize = '24')
plt.legend(fontsize = 'x-large', loc = 'best')
plt.savefig("{}_rdf.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_rdf.svg".format(args.o), format='svg', dpi=300)     
