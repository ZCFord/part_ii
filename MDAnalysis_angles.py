import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import math
import argparse

#rcParams['axes.labelsize'] = 8
#rcParams['xtick.labelsize'] = 8
#rcParams['ytick.labelsize'] = 8
#rcParams['legend.fontsize'] = 10
#rcParams['font.family'] = ['sans-serif']
#rcParams['font.sans-serif'] = ['Arial']
rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'
#matplotlib.rcParams['axes.unicode_minus'] = False
#rcParams['figure.subplot.wspace']= 0.3
#rcParams['figure.subplot.bottom']= 0.1
#rcParams['figure.subplot.hspace']= 0.2
#rcParams['figure.subplot.left']= 0.125
#rcParams['figure.subplot.right']= 0.9
#rcParams['figure.subplot.top']= 0.9
#matplotlib.rcParams['axes.linewidth']= .5

##########
# ARGUMENT HANDLING
##########

parser = argparse.ArgumentParser('Calculates angle of M2, M3, and M4 helices')

parser.add_argument(
	"-f",
	type = str,
	default = './mdord.pdb',
	help="Structure (pdb)"
)

parser.add_argument(
	"-t",
	type = str,	
	default = 'md-c.xtc',
	help="Trajectory (xtc)"
)

parser.add_argument(
	"-o",
	type = str,
	default = None,
	help="output name"
)

args = parser.parse_args()

##########
# DEFINITIONS
##########

def angle (sel1, sel2, sel3) : #input three residues
	res_i = protein.select_atoms(sel1 + ' and name CA').positions
	res_j = protein.select_atoms(sel2 + ' and name CA').positions
	res_k = protein.select_atoms(sel3 + ' and name CA').positions
	angle1 = res_i - res_j
	angle2 = res_k - res_j
	dot = np.dot(angle1[0], angle2[0]) #dot product
	cosTH = dot / (np.linalg.norm(angle1) * np.linalg.norm(angle2))
	return  math.degrees(np.arccos(cosTH))

def dis (sel1, sel2, x) :
	res_i = x.select_atoms(sel1 + ' and name CA').positions
	res_j = x.select_atoms(sel2 + ' and name CA').positions
	dis = np.linalg.norm(res_i - res_j)
	return dis

##########
# CALCULATIONS
##########

u = MDAnalysis.Universe(args.f, args.t)
protein = u.select_atoms('protein')
m2_A = []
m2_B = []
m3_A = []
m3_B = []
m4_A = []
m4_B = []

for ts in u.trajectory:
	m2_A.append(angle('resid 185 and segid A', 'resid 198 and segid A', 'resid 227 and segid A'))
	m2_B.append(angle('resid 185 and segid B', 'resid 198 and segid B', 'resid 227 and segid B'))
	m3_A.append(angle('resid 262 and segid A', 'resid 254 and segid A', 'resid 233 and segid A'))
	m3_B.append(angle('resid 262 and segid B', 'resid 254 and segid B', 'resid 233 and segid B'))
	m4_A.append(angle('resid 308 and segid A', 'resid 322 and segid A', 'resid 332 and segid A'))
	m4_B.append(angle('resid 308 and segid B', 'resid 322 and segid B', 'resid 332 and segid B'))

M1_A = []
M1_B = []
M2_A = []
M2_B = []
M3_A = []
M3_B = []
M4_A = []
M4_B = []

for ts in u.trajectory:
	M1_A.append(dis('resid 97 and segid A', 'resid 75 and segid A', protein))
	M1_B.append(dis('resid 97 and segid B', 'resid 75 and segid B', protein))
	M2_A.append(dis('resid 195 and segid A', 'resid 222 and segid A', protein))
	M2_B.append(dis('resid 195 and segid B', 'resid 222 and segid B', protein))
	M3_A.append(dis('resid 262 and segid A', 'resid 236 and segid A', protein))
	M3_B.append(dis('resid 262 and segid B', 'resid 236 and segid B', protein))
	M4_A.append(dis('resid 302 and segid A', 'resid 328 and segid A', protein))
	M4_B.append(dis('resid 302 and segid B', 'resid 328 and segid B', protein))

M1_a = []
M1_b = []
M2_a = []
M2_b = []
M3_a = []
M3_b = []
M4_a = []
M4_b = []

for ts in u.trajectory:
	M1_a.append((protein.select_atoms('resid 97 and segid A and name CA').positions)[0][2] - (protein.select_atoms('resid 75 and segid A and name CA').positions)[0][2])
	M1_b.append((protein.select_atoms('resid 97 and segid B and name CA').positions)[0][2] - (protein.select_atoms('resid 75 and segid B and name CA').positions)[0][2])
	M2_a.append((protein.select_atoms('resid 195 and segid A and name CA').positions)[0][2] - (protein.select_atoms('resid 222 and segid A and name CA').positions)[0][2])
	M2_b.append((protein.select_atoms('resid 195 and segid B and name CA').positions)[0][2] - (protein.select_atoms('resid 222 and segid B and name CA').positions)[0][2])
	M3_a.append((protein.select_atoms('resid 262 and segid A and name CA').positions)[0][2] - (protein.select_atoms('resid 236 and segid A and name CA').positions)[0][2])
	M3_b.append((protein.select_atoms('resid 262 and segid B and name CA').positions)[0][2] - (protein.select_atoms('resid 236 and segid B and name CA').positions)[0][2])
	M4_a.append((protein.select_atoms('resid 302 and segid A and name CA').positions)[0][2] - (protein.select_atoms('resid 328 and segid A and name CA').positions)[0][2])
	M4_b.append((protein.select_atoms('resid 302 and segid B and name CA').positions)[0][2] - (protein.select_atoms('resid 328 and segid B and name CA').positions)[0][2])

M1_A_angle = []
M1_B_angle = []
M2_A_angle = []
M2_B_angle = []
M3_A_angle = []
M3_B_angle = []
M4_A_angle = []
M4_B_angle = []

for i in M1_A:
   	cosTH = M1_a[(M1_A.index(i))-1]/i
	M1_A_angle.append(math.degrees(np.arccos(cosTH)))
for i in M1_B:
   	cosTH = M1_b[(M1_B.index(i))-1]/i
	M1_B_angle.append(math.degrees(np.arccos(cosTH)))
for i in M2_A:
   	cosTH = M2_a[(M2_A.index(i))-1]/i
	M2_A_angle.append(math.degrees(np.arccos(cosTH)))
for i in M2_B:
   	cosTH = M2_b[(M2_B.index(i))-1]/i
	M2_B_angle.append(math.degrees(np.arccos(cosTH)))
for i in M3_A:
   	cosTH = M3_a[(M3_A.index(i))-1]/i
	M3_A_angle.append(math.degrees(np.arccos(cosTH)))
for i in M3_B:
   	cosTH = M3_b[(M3_B.index(i))-1]/i
	M3_B_angle.append(math.degrees(np.arccos(cosTH)))
for i in M4_A:
   	cosTH = M4_a[(M4_A.index(i))-1]/i
	M4_A_angle.append(math.degrees(np.arccos(cosTH)))
for i in M4_B:
   	cosTH = M4_b[(M4_B.index(i))-1]/i
	M4_B_angle.append(math.degrees(np.arccos(cosTH)))
	
##########
# PLOTTING
##########

plt.figure(figsize=(30,10))
A = plt.subplot(131)
plt.plot(m2_A, label = 'Chain A', color = 'blue', linewidth = 2)
plt.plot(m2_B, label = 'Chain B', color = 'gray', linewidth = 2)
plt.setp(A.get_xticklabels(), fontsize = 18)
plt.setp(A.get_yticklabels(), fontsize = 18)
plt.title('M2 Helix', fontsize = 36)
v = [0,100,0,360]
plt.axis(v)
plt.xlabel("time (ns)", fontsize = '24')
plt.ylabel("Angle (degrees)", fontsize = '24')
B = plt.subplot(132, sharex = A, sharey = A)
plt.plot(m3_A, label = 'Chain A', color = 'orange', linewidth = 2)
plt.plot(m3_B, label = 'Chain B', color = 'gray', linewidth = 2)
plt.setp(B.get_xticklabels(), fontsize = 18)
plt.setp(B.get_yticklabels(), fontsize = 18)
plt.xlabel("time (ns)", fontsize = '24')
plt.title('M3 Helix', fontsize = 36)                                         
C = plt.subplot(133, sharex = A, sharey = A)
plt.plot(m4_A, label = 'Chain A', color = 'green', linewidth = 2)
plt.plot(m4_B, label = 'Chain B', color = 'gray', linewidth = 2)
plt.setp(C.get_xticklabels(), fontsize = 18)
plt.setp(C.get_yticklabels(), fontsize = 18)
plt.xlabel("time (ns)", fontsize = '24')
plt.title('M4 Helix', fontsize = 36)
plt.legend(fontsize = 'xx-large', loc = 'best')  
plt.savefig("{}_angles.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_angles.svg".format(args.o), format='svg', dpi=300)

plt.clf()
plt.figure(figsize=(20,20))
one = plt.subplot(221)
plt.plot(M1_A_angle, label = 'Chain A', color = 'purple', linewidth = 2, alpha = 0.5)
plt.plot(M1_B_angle, label = 'Chain B', color = 'purple', linewidth = 2, linestyle = '--')
plt.setp(one.get_xticklabels(), fontsize = 14)
plt.setp(one.get_yticklabels(), fontsize = 14)
plt.xlabel("time (ns)", fontsize = '18')
plt.ylabel("Tilt Angle (degrees)", fontsize = '18')
plt.title('M1 Helix', fontsize = 24)
plt.legend(fontsize = 'x-large', loc = 'best')
v = [0,100,0,80]
plt.axis(v)
two = plt.subplot(222, sharex = one, sharey = one)
plt.plot(M2_A_angle, label = 'Chain A', color = 'blue', linewidth = 2, alpha = 0.5)
plt.plot(M2_B_angle, label = 'Chain B', color = 'blue', linewidth = 2, linestyle = '--')
plt.setp(two.get_xticklabels(), fontsize = 14)
plt.setp(two.get_yticklabels(), fontsize = 14)
plt.xlabel("time (ns)", fontsize = '18')
plt.ylabel("Tilt Angle (degrees)", fontsize = '18')
plt.title('M2 Helix', fontsize = 24)
three = plt.subplot(223, sharex = one, sharey = one)
plt.plot(M3_A_angle, label = 'Chain A', color = 'orange', linewidth = 2, alpha = 0.5)
plt.plot(M3_B_angle, label = 'Chain B', color = 'orange', linewidth = 2, linestyle = '--')
plt.setp(three.get_xticklabels(), fontsize = 14)
plt.setp(three.get_yticklabels(), fontsize = 14)
plt.xlabel("time (ns)", fontsize = '18')
plt.ylabel("Tilt Angle (degrees)", fontsize = '18')
plt.title('M3 Helix', fontsize = 24)
four = plt.subplot(224, sharex = one, sharey = one)
plt.plot(M4_A_angle, label = 'Chain A', color = 'green', linewidth = 2, alpha = 0.5)
plt.plot(M4_B_angle, label = 'Chain B', color = 'green', linewidth = 2, linestyle = '--')
plt.setp(four.get_xticklabels(), fontsize = 14)
plt.setp(four.get_yticklabels(), fontsize = 14)
plt.xlabel("time (ns)", fontsize = '18')
plt.ylabel("Tilt Angle (degrees)", fontsize = '18')
plt.title('M4 Helix', fontsize = 24)
plt.suptitle('Tilt Angle relative to z axis', fontsize = 36)
plt.savefig("{}_newangles.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_newangles.svg".format(args.o), format='svg', dpi=300)
