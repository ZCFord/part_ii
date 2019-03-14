import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False

DCPC = np.loadtxt('./DCPC/md/most_frames/stress/stress.txt')
DOPC = np.loadtxt('./DOPC/md/most_frames/stress/stress.txt')
DPPC = np.loadtxt('./DPPC/md/most_frames/stress/stress.txt')
DWPC = np.loadtxt('./DWPC/md/most_frames/stress/stress.txt')
POPC = np.loadtxt('./POPC/md/most_frames/stress/stress.txt')

c = DCPC[:,0]-(DCPC[-1,0]/2)
o = DOPC[:,0]-(DOPC[-1,0]/2)
p = DPPC[:,0]-(DPPC[-1,0]/2)
w = DWPC[:,0]-(DWPC[-1,0]/2)
pop = POPC[:,0]-(POPC[-1,0]/2)

plt.figure(figsize=(10,8))
plt.title('Bilayer Pressure Profile', fontsize = 30)
plt.plot(c*10, (DCPC[:,1]+DCPC[:,5])/2, label = '3-bead Tail (CDC)', color = 'blue')
plt.plot(o*10, (DOPC[:,1]+DOPC[:,5])/2, label = '4-bead Tail (CDCC)', color = 'purple')
plt.plot(p*10, (DPPC[:,1]+DPPC[:,5])/2, label = '4-bead Tail (CCCC)', color = 'red')
plt.plot(pop*10, (POPC[:,1]+POPC[:,5])/2, label = '4-bead Tail (CDCC/CCCC)', color = 'orange')
plt.plot(w*10, (DWPC[:,1]+DWPC[:,5])/2, label = '5-bead Tail (CDCCC)', color = 'green')
v = [-50,50,-500,600]
plt.axis(v)
plt.ylabel("P (bar)", fontsize = '24')
plt.xlabel("Z ($\AA$)", fontsize = '24')
plt.legend(fontsize = 'x-large', loc = 'best')
plt.savefig("stress_all.png", format='png', dpi=300)
plt.savefig("stress_all.svg", format='svg', dpi=300)   

