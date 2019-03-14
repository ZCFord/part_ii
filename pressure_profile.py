import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False

file = np.loadtxt('./stress.txt')

plt.figure(figsize=(10,8))
plt.title('Bilayer Pressure Profile', fontsize = 30)
plt.plot(file[:,0]*10, (file[:,1]+file[:,5])/2, label = '$P_L = -(\sigma_{xx} + \sigma_{yy})/2$')
plt.plot(file[:,0]*10, file[:,9], label = '$P_N = -\sigma_{zz}$')
v = [0,133,-500,600]
plt.axis(v)
plt.ylabel("P (bar)", fontsize = '24')
plt.xlabel("Z ($\AA$)", fontsize = '24')
plt.legend(fontsize = 'x-large', loc = 'best')
plt.savefig("stress.png", format='png', dpi=300)
plt.savefig("stress.svg", format='svg', dpi=300)   

