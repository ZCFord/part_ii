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

parser = argparse.ArgumentParser('Plot APL')

parser.add_argument(
	"-o",
	type = str,
	default = None,
	help = "output name"
)

args = parser.parse_args()

APL = np.loadtxt('./apl.xvg', skiprows = 15)

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

av = running_mean(APL[:,1]*100, 10)
up = running_mean(APL[:,2]*100, 10)
lo = running_mean(APL[:,3]*100, 10)

plt.figure(figsize=(10,8))
x = np.arange(4,96,1)
plt.plot(APL[:,0]/1000, APL[:,1]*100, label = 'Average APL', color = 'blue', linewidth = 2)
plt.plot(x, av, label = '10ns Running Average', color = 'blue', linewidth = 1, linestyle = 'dashed')
plt.plot(x, up, label = 'Upper Leaflet', color = 'red', linewidth = 1)
plt.plot(x, lo, label = 'Lower Leaflet', color = 'green', linewidth = 1)
v = [0,100,50,120]
plt.axis(v)
plt.ylabel("Area per Lipid ($\AA^2$)", fontsize = '24')
plt.xlabel("Time (ns)", fontsize = '24')
plt.legend(fontsize = 'x-large', loc = 'best')
plt.savefig("{}_apl.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_apl.svg".format(args.o), format='svg', dpi=300)
