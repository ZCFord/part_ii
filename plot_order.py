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

parser = argparse.ArgumentParser('Plots the order parameter of the lipids')

parser.add_argument(
	"-n",
	type = float,
	required = True,
	help = "Number of carbons in lipid tail"
)

parser.add_argument(
	"-d",
	default = "Yes",
	help = "Is there a double bond?"
)

parser.add_argument(
	"-o",
	type = str,
	default = None,
	help = "output name"
)

args = parser.parse_args()

U1 = np.loadtxt('./upper_graph_1.xvg', skiprows = 17)
U2 = np.loadtxt('./upper_graph_2.xvg', skiprows = 17)
L1 = np.loadtxt('./lower_graph_1.xvg', skiprows = 17)
L2 = np.loadtxt('./lower_graph_2.xvg', skiprows = 17)

plt.figure(figsize=(24,8))
UP = plt.subplot(121)
plt.plot(U1[:,0]+1, U1[:,1], label = 'sn1 Chain', color = 'blue')
plt.plot(U2[:,0]+1, U2[:,1], label = 'sn2 Chain', color = 'red')
plt.xlabel('Carbon Atom', fontsize = '18')
plt.ylabel('Deuterium Order Parameter (Scd)', fontsize = '18')
plt.setp(UP.get_xticklabels(), fontsize = 18)
plt.setp(UP.get_yticklabels(), fontsize = 18)
plt.title("Upper Leaflet", fontsize = 36)
if args.d == "Yes" :
	plt.axvspan(9, 10, facecolor = '0.5', alpha = 0.2)
LO = plt.subplot(122, sharex = UP, sharey = UP)
plt.plot(L1[:,0]+1, L1[:,1], label = 'sn1 Chain', color = 'blue')
plt.plot(L2[:,0]+1, L2[:,1], label = 'sn2 Chain', color = 'red')
plt.setp(LO.get_xticklabels(), fontsize = 18)
plt.setp(LO.get_yticklabels(), fontsize = 18)
plt.xlabel('Carbon Atom', fontsize = '18')
plt.ylabel('Deuterium Order Parameter (Scd)', fontsize = '18')
plt.legend(fontsize = 'x-large', loc = 'best')
plt.title("Lower Leaflet", fontsize = 36)
if args.d == "Yes" :
	plt.axvspan(9, 10, facecolor = '0.5', alpha = 0.2)
v = [1, args.n, 0, 0.25]
plt.axis(v)
plt.savefig("{}_orderparam.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_orderparam.svg".format(args.o), format='svg', dpi=300)
