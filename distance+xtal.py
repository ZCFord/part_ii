import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
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

parser = argparse.ArgumentParser('Calculates FEZ distances')

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

def dis (sel1, sel2, x) :
	res_i = x.select_atoms(sel1 + ' and name CA').coordinates()
	res_j = x.select_atoms(sel2 + ' and name CA').coordinates()
	dis = np.linalg.norm(res_i - res_j)
	return dis

##########
# CALCULATIONS
##########

u = MDAnalysis.Universe(args.f, args.t)
protein = u.select_atoms('protein')

up = MDAnalysis.Universe('./UP_correctpro.pdb')
up_pro = up.select_atoms('protein')

dn = MDAnalysis.Universe('./DN_correctH.pdb')
dn_pro = dn.select_atoms('protein')

Fen_A = []
Fen_B = []
Zip_A = []
Zip_B = []
Exp_A = []
Exp_B = []

for ts in u.trajectory:
	Fen_A.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein))
	Fen_B.append(dis('resid 324 and segid B', 'resid 198 and segid A', protein))
	Zip_A.append(dis('resid 326 and segid A', 'resid 237 and segid A', protein))
	Zip_B.append(dis('resid 326 and segid B', 'resid 237 and segid B', protein))
	Exp_A.append(dis('resid 322 and segid A', 'resid 212 and segid A', protein))
	Exp_B.append(dis('resid 322 and segid B', 'resid 212 and segid B', protein))	

FAup = []
FBup = []
ZAup = []
ZBup = []
EAup = []
EBup = []

FAup.append(dis('resid 324 and segid A', 'resid 198 and segid B', up_pro))
FBup.append(dis('resid 324 and segid B', 'resid 198 and segid A', up_pro))
ZAup.append(dis('resid 326 and segid A', 'resid 237 and segid A', up_pro))
ZBup.append(dis('resid 326 and segid B', 'resid 237 and segid B', up_pro))
EAup.append(dis('resid 322 and segid A', 'resid 212 and segid A', up_pro))
EBup.append(dis('resid 322 and segid B', 'resid 212 and segid B', up_pro))

FAdn = []
FBdn = []
ZAdn = []
ZBdn = []
EAdn = []
EBdn = []

FAdn.append(dis('resid 324 and segid A', 'resid 198 and segid B', dn_pro))
FBdn.append(dis('resid 324 and segid B', 'resid 198 and segid A', dn_pro))
ZAdn.append(dis('resid 326 and segid A', 'resid 237 and segid A', dn_pro))
ZBdn.append(dis('resid 326 and segid B', 'resid 237 and segid B', dn_pro))
EAdn.append(dis('resid 322 and segid A', 'resid 212 and segid A', dn_pro))
EBdn.append(dis('resid 322 and segid B', 'resid 212 and segid B', dn_pro))

##########
# PLOTTING
##########

#plt.figure(figsize=(20,16))
#A = plt.subplot(211)
#plt.plot(Fen_A, label = 'Fenestration', color = 'blue', linewidth = 3)
#plt.plot(Zip_A, label = 'Zipper', color = 'green', linewidth = 3)
#plt.plot(Exp_A, label = 'Expansion', color = 'red', linewidth = 3)
#plt.setp(A.get_xticklabels(), visible=False)
#plt.setp(A.get_yticklabels(), fontsize = 18)
#plt.ylabel("z ($\AA$)", fontsize = '24')
#plt.legend(fontsize = 'xx-large', loc = 'best')    
#plt.title('Chain A', fontsize = 36)                                     
#B = plt.subplot(212, sharex = A, sharey = A)
#plt.plot(Fen_B, label = 'Fenestration', color = 'blue', linewidth = 3)
#plt.plot(Zip_B, label = 'Zipper', color = 'green', linewidth = 3)
#plt.plot(Exp_B, label = 'Expansion', color = 'red', linewidth = 3)
#plt.setp(B.get_xticklabels(), fontsize = 18)
#plt.setp(B.get_yticklabels(), fontsize = 18)
#plt.xlabel("time (ns)", fontsize = '24')
#plt.ylabel("z ($\AA$)", fontsize = '24')
#v = [0,100,3,18]
#plt.axis(v)
#plt.title('Chain B', fontsize = 36) 
#plt.savefig("{}_FEZ_AB.png".format(args.o), format='png', dpi=300)
#plt.savefig("{}_FEZ_AB.svg".format(args.o), format='svg', dpi=300)

plt.clf()
plt.figure(figsize=(16,8))
FZ = plt.subplot(121)
plt.scatter(Zip_A, Fen_A, label = 'Chain A', color = 'green')
plt.scatter(Zip_B, Fen_B, label = 'Chain B', color = 'gray')
plt.plot(Zip_A, Fen_A, label = 'Chain A', color = 'green', alpha = 0.5)
plt.plot(Zip_B, Fen_B, label = 'Chain B', color = 'gray', alpha = 0.5)
plt.scatter(ZAup, FAup, label = '4BW5 (Up)', color = 'black', marker = 'D', s = 40)
plt.scatter(ZAdn, FAdn, label = '4XDJ (Down)', color = 'blue', marker = 'D', s = 40)
plt.xlabel('Zipper ($\AA$)', fontsize = '18')
plt.ylabel('Fenestration ($\AA$)', fontsize = '18')
plt.setp(FZ.get_xticklabels(), fontsize = 18)
plt.setp(FZ.get_yticklabels(), fontsize = 18)
plt.legend(fontsize = 'x-large', loc = 'best')
FE = plt.subplot(122, sharex = FZ, sharey = FZ)
plt.scatter(Exp_A, Fen_A, label = 'Chain A', color = 'red')
plt.scatter(Exp_B, Fen_B, label = 'Chain B', color = 'pink')
plt.plot(Exp_A, Fen_A, label = 'Chain A', color = 'red', alpha = 0.5)
plt.plot(Exp_B, Fen_B, label = 'Chain B', color = 'pink', alpha = 0.5)
plt.scatter(EAup, FAup, label = '4BW5 (Up)', color = 'black', marker = 'D', s = 40)
plt.scatter(EAdn, FAdn, label = '4XDJ (Down)', color = 'blue', marker = 'D', s = 40)
plt.xlabel('Expansion ($\AA$)', fontsize = '18')
plt.setp(FE.get_xticklabels(), fontsize = 18)
plt.setp(FE.get_yticklabels(), visible = False)
v = [3,18,3,16]
plt.axis(v)
plt.legend(fontsize = 'x-large', loc = 'best')
plt.savefig("{}_FZ_FE.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_FZ_FE.svg".format(args.o), format='svg', dpi=300)
