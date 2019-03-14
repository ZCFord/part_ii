import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import argparse

rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'

parser = argparse.ArgumentParser('Calculates angle of M2, M3, and M4 helices')

parser.add_argument(
	"-o",
	type = str,
	default = None,
	help="output name"
)

args = parser.parse_args()

A2 = np.loadtxt('./2a.xvg', skiprows = 17)
B2 = np.loadtxt('./2b.xvg', skiprows = 17)
A3 = np.loadtxt('./3a.xvg', skiprows = 17)
B3 = np.loadtxt('./3b.xvg', skiprows = 17)
A4 = np.loadtxt('./4a.xvg', skiprows = 17)
B4 = np.loadtxt('./4b.xvg', skiprows = 17)

a2 = []

for i in A2[:,1]:
    a2.append(i - A2[0,1])
    
b2 = []

for i in B2[:,1]:
    b2.append(i - B2[0,1])
    
a3 = []

for i in A3[:,1]:
    a3.append(i - A3[0,1])
    
b3 = []

for i in B3[:,1]:
    b3.append(i - B3[0,1])
    
a4 = []

for i in A4[:,1]:
    a4.append(i - A4[0,1])
    
b4 = []

for i in B4[:,1]:
    b4.append(i - B4[0,1])

plt.figure(figsize=(30,10))
A = plt.subplot(131)
plt.plot(A2[:,0], a2, label = 'Chain A', color = 'blue', linewidth = 2)
plt.plot(B2[:,0], b2, label = 'Chain B', color = 'gray', linewidth = 2)
plt.setp(A.get_xticklabels(), fontsize = 18)
plt.setp(A.get_yticklabels(), fontsize = 18)
plt.title('M2 Helix', fontsize = 30)
v = [0,100,-10,10]
plt.axis(v)
plt.xlabel("time (ns)", fontsize = '24')
plt.ylabel("Angle (degrees)", fontsize = '24')
B = plt.subplot(132, sharex = A, sharey = A)
plt.plot(A3[:,0], a3, label = 'Chain A', color = 'orange', linewidth = 2)
plt.plot(B3[:,0], b2, label = 'Chain B', color = 'gray', linewidth = 2)
plt.setp(B.get_xticklabels(), fontsize = 18)
plt.setp(B.get_yticklabels(), fontsize = 18)
plt.xlabel("time (ns)", fontsize = '24')
plt.title('M3 Helix', fontsize = 30)                                         
C = plt.subplot(133, sharex = A, sharey = A)
plt.plot(A4[:,0], a4, label = 'Chain A', color = 'green', linewidth = 2)
plt.plot(B4[:,0], b4, label = 'Chain B', color = 'gray', linewidth = 2)
plt.setp(C.get_xticklabels(), fontsize = 18)
plt.setp(C.get_yticklabels(), fontsize = 18)
plt.xlabel("time (ns)", fontsize = '24')
plt.title('M4 Helix', fontsize = 30)
plt.legend(fontsize = 'xx-large', loc = 'best')
plt.suptitle('Change in helix angle over time', fontsize = 36)  
plt.savefig("{}_lol-angles.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_lol-angles.svg".format(args.o), format='svg', dpi=300)
