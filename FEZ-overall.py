import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from matplotlib import rc, rcParams

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False

##########
# DEFINITIONS
##########

def dis (sel1, sel2, x) :
	res_i = x.select_atoms(sel1 + ' and name CA').positions
	res_j = x.select_atoms(sel2 + ' and name CA').positions
	dis = np.linalg.norm(res_i - res_j)
	return dis

def FEZ_traj(A, B, C, D) :
	u = MDAnalysis.Universe(A, B)
	Fen_A = []
	Fen_B = []
	Zip_A = []
	Zip_B = []
	Exp_A = []
	Exp_B = []
	Fen_A_2 = []
	Fen_B_2 = []
	Zip_A_2 = []
	Zip_B_2 = []
	Exp_A_2 = []
	Exp_B_2 = []
	for ts in u.trajectory:
		protein = u.select_atoms('protein')
		Fen_A.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein))
		Fen_B.append(dis('resid 324 and segid B', 'resid 198 and segid A', protein))
		Zip_A.append(dis('resid 326 and segid A', 'resid 237 and segid A', protein))
		Zip_B.append(dis('resid 326 and segid B', 'resid 237 and segid B', protein))
		Exp_A.append(dis('resid 322 and segid A', 'resid 212 and segid A', protein))
		Exp_B.append(dis('resid 322 and segid B', 'resid 212 and segid B', protein))	
	for i in Fen_A :
		if i < 5.5 :
			Fen_A_2.append(Fen_A.index(i))
	for i in Fen_B :
		if i < 5.5 :
			Fen_B_2.append(Fen_B.index(i))
	for i in Zip_A :
		if i > 11.5 :
			Zip_A_2.append(Zip_A.index(i))
	for i in Zip_B :
		if i > 11.5 :
			Zip_B_2.append(Zip_B.index(i))
	for i in Exp_A :
		if i > 11.5 :
			Exp_A_2.append(Exp_A.index(i))
	for i in Exp_B :
		if i > 11.5 :
			Exp_B_2.append(Exp_B.index(i))
	for i in Fen_A_2 :
		if i in Zip_A_2 and i in Exp_A_2 :
			C.append(i)
	for i in Fen_B_2 :
		if i in Zip_B_2 and i in Exp_B_2 :
			D.append(i)
	if C == [] :
		C.append(float(101))
	if D == [] :
		D.append(float(101))

##########
# CALCULATIONS
##########

##### DFPC
### 1b

F1b1chA = []
F1b1chB = []
FEZ_traj('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc', F1b1chA, F1b1chB)

F1b2chA = []
F1b2chB = []
FEZ_traj('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', F1b2chA, F1b2chB)

### -30b

F30b1chA = []
F30b1chB = []
FEZ_traj('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', F30b1chA, F30b1chB)

F30b2chA = []
F30b2chB = []
FEZ_traj('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', F30b2chA, F30b2chB)

### -40b

F40b1chA = []
F40b1chB = []
FEZ_traj('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', F40b1chA, F40b1chB)

F40b2chA = []
F40b2chB = []
FEZ_traj('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', F40b2chA, F40b2chB)

F40b3chA = []
F40b3chB = []
FEZ_traj('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', F40b3chA, F40b3chB)

##### DOPC
### 1b

O1b1chA = []
O1b1chB = []
FEZ_traj('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc', O1b1chA, O1b1chB)

O1b2chA = []
O1b2chB = []
FEZ_traj('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', O1b2chA, O1b2chB)

### -30b

O30b1chA = []
O30b1chB = []
FEZ_traj('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', O30b1chA, O30b1chB)

O30b2chA = []
O30b2chB = []
FEZ_traj('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', O30b2chA, O30b2chB)

### -40b

O40b1chA = []
O40b1chB = []
FEZ_traj('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', O40b1chA, O40b1chB)

O40b2chA = []
O40b2chB = []
FEZ_traj('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', O40b2chA, O40b2chB)

O40b3chA = []
O40b3chB = []
FEZ_traj('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', O40b3chA, O40b3chB)

### -50b

O50b1chA = []
O50b1chB = []
FEZ_traj('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', O50b1chA, O50b1chB)

O50b2chA = []
O50b2chB = []
FEZ_traj('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', O50b2chA, O50b2chB)

##### DVPC
### 1b

V1b1chA = []
V1b1chB = []
FEZ_traj('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc', V1b1chA, V1b1chB)

V1b2chA = []
V1b2chB = []
FEZ_traj('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', V1b2chA, V1b2chB)

### -40b

V40b1chA = []
V40b1chB = []
FEZ_traj('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', V40b1chA, V40b1chB)

V40b2chA = []
V40b2chB = []
FEZ_traj('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', V40b2chA, V40b2chB)

V40b3chA = []
V40b3chB = []
FEZ_traj('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', V40b3chA, V40b3chB)

### -50b

V50b1chA = []
V50b1chB = []
FEZ_traj('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', V50b1chA, V50b1chB)

V50b2chA = []
V50b2chB = []
FEZ_traj('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', V50b2chA, V50b2chB)

### -60b

V60b1chA = []
V60b1chB = []
FEZ_traj('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', V60b1chA, V60b1chB)

V60b2chA = []
V60b2chB = []
FEZ_traj('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', V60b2chA, V60b2chB)

##### POPC
### 1b

P1b1chA = []
P1b1chB = []
FEZ_traj('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc', P1b1chA, P1b1chB)

P1b2chA = []
P1b2chB = []
FEZ_traj('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc', P1b2chA, P1b2chB)

### -50b

P50b1chA = []
P50b1chB = []
FEZ_traj('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', P50b1chA, P50b1chB)

P50b2chA = []
P50b2chB = []
FEZ_traj('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', P50b2chA, P50b2chB)

##########
# PLOTTING
##########

fig = plt.figure(figsize=(16,6))
ax = plt.subplot(111)

### Chain A

plt.axvline(F1b1chA[0], 0, 1, linewidth = 2, color = 'blue', label = 'DFPC 1b Run 1 (Chain A)')
plt.axvline(F1b2chA[0], 0, 1, linewidth = 2, color = 'blue', alpha = 0.66, label = 'DFPC 1b Run 2 (Chain A)')
plt.axvline(F30b1chA[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '--', label = 'DFPC 30b Run 1 (Chain A)')
plt.axvline(F30b2chA[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '--', alpha = 0.66, label = 'DFPC 30b Run 2 (Chain A)')
plt.axvline(F40b1chA[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '-.', label = 'DFPC 40b Run 1 (Chain A)')
plt.axvline(F40b2chA[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '-.', alpha = 0.66, label = 'DFPC 40b Run 2 (Chain A)')
plt.axvline(F40b3chA[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '-.', alpha = 0.33, label = 'DFPC 40b Run 3 (Chain A)')
plt.axvline(O1b1chA[0], 0, 1, linewidth = 2, color = 'green', label = 'DOPC 1b Run 1 (Chain A)')
plt.axvline(O1b2chA[0], 0, 1, linewidth = 2, color = 'green', alpha = 0.66, label = 'DOPC 1b Run 2 (Chain A)')
plt.axvline(O30b1chA[0], 0, 1, linewidth = 2, color = 'green', linestyle = '--', label = 'DOPC 30b Run 1 (Chain A)')
plt.axvline(O30b2chA[0], 0, 1, linewidth = 2, color = 'green', linestyle = '--', alpha = 0.66, label = 'DOPC 30b Run 2 (Chain A)')
plt.axvline(O40b1chA[0], 0, 1, linewidth = 2, color = 'green', linestyle = '-.', label = 'DOPC 40b Run 1 (Chain A)')
plt.axvline(O40b2chA[0], 0, 1, linewidth = 2, color = 'green', linestyle = '-.', alpha = 0.66, label = 'DOPC 40b Run 2 (Chain A)')
plt.axvline(O40b3chA[0], 0, 1, linewidth = 2, color = 'green', linestyle = '-.', alpha = 0.33, label = 'DOPC 1b Run 3 (Chain A)')
plt.axvline(O50b1chA[0], 0, 1, linewidth = 2, color = 'green', linestyle = ':', label = 'DOPC 50b Run 1 (Chain A)')
plt.axvline(O50b2chA[0], 0, 1, linewidth = 2, color = 'green', linestyle = ':', alpha = 0.66, label = 'DOPC 50b Run 2 (Chain A)')
plt.axvline(V1b1chA[0], 0, 1, linewidth = 2, color = 'red', label = 'DVPC 1b Run 1 (Chain A)')
plt.axvline(V1b2chA[0], 0, 1, linewidth = 2, color = 'red', alpha = 0.66, label = 'DVPC 1b Run 2 (Chain A)')
plt.axvline(V40b1chA[0], 0, 1, linewidth = 2, color = 'red', linestyle = '--', label = 'DVPC 40b Run 1 (Chain A)')
plt.axvline(V40b2chA[0], 0, 1, linewidth = 2, color = 'red', linestyle = '--', alpha = 0.66, label = 'DVPC 40b Run 2 (Chain A)')
plt.axvline(V40b3chA[0], 0, 1, linewidth = 2, color = 'red', linestyle = '--', alpha = 0.33, label = 'DVPC 40b Run 3 (Chain A)')
plt.axvline(V50b1chA[0], 0, 1, linewidth = 2, color = 'red', linestyle = '-.', label = 'DVPC 50b Run 1 (Chain A)')
plt.axvline(V50b2chA[0], 0, 1, linewidth = 2, color = 'red', linestyle = '-.', alpha = 0.66, label = 'DVPC 50b Run 2 (Chain A)')
plt.axvline(V60b1chA[0], 0, 1, linewidth = 2, color = 'red', linestyle = ':', label = 'DVPC 60b Run 1 (Chain A)')
plt.axvline(V60b2chA[0], 0, 1, linewidth = 2, color = 'red', linestyle = ':', alpha = 0.66, label = 'DVPC 60b Run 2 (Chain A)')
plt.axvline(P1b1chA[0], 0, 1, linewidth = 2, color = 'black', label = 'POPC 1b Run 1 (Chain A)')
plt.axvline(P1b2chA[0], 0, 1, linewidth = 2, color = 'black', alpha = 0.66, label = 'POPC 1b Run 2 (Chain A)')
plt.axvline(P50b1chA[0], 0, 1, linewidth = 2, color = 'black', linestyle = '--', label = 'POPC 50b Run 1 (Chain A)')
plt.axvline(P50b2chA[0], 0, 1, linewidth = 2, color = 'black', linestyle = '--', alpha = 0.66, label = 'POPC 50b Run 2 (Chain A)')

### Chain B

plt.axvline(F1b1chB[0], 0, 1, linewidth = 2, color = 'cyan', label = 'DFPC 1b Run 1 (Chain B)')
plt.axvline(F1b2chB[0], 0, 1, linewidth = 2, color = 'cyan', alpha = 0.66, label = 'DFPC 1b Run 2 (Chain B)')
plt.axvline(F30b1chB[0], 0, 1, linewidth = 2, color = 'cyan', linestyle = '--', label = 'DFPC 30b Run 1 (Chain B)')
plt.axvline(F30b2chB[0], 0, 1, linewidth = 2, color = 'cyan', linestyle = '--', alpha = 0.66, label = 'DFPC 30b Run 2 (Chain B)')
plt.axvline(F40b1chB[0], 0, 1, linewidth = 2, color = 'cyan', linestyle = '-.', label = 'DFPC 40b Run 1 (Chain B)')
plt.axvline(F40b2chB[0], 0, 1, linewidth = 2, color = 'cyan', linestyle = '-.', alpha = 0.66, label = 'DFPC 40b Run 2 (Chain B)')
plt.axvline(F40b3chB[0], 0, 1, linewidth = 2, color = 'cyan', linestyle = '-.', alpha = 0.33, label = 'DFPC 40b Run 3 (Chain B)')
plt.axvline(O1b1chB[0], 0, 1, linewidth = 2, color = 'yellow', label = 'DOPC 1b Run 1 (Chain B)')
plt.axvline(O1b2chB[0], 0, 1, linewidth = 2, color = 'yellow', alpha = 0.66, label = 'DOPC 1b Run 2 (Chain B)')
plt.axvline(O30b1chB[0], 0, 1, linewidth = 2, color = 'yellow', linestyle = '--', label = 'DOPC 30b Run 1 (Chain B)')
plt.axvline(O30b2chB[0], 0, 1, linewidth = 2, color = 'yellow', linestyle = '--', alpha = 0.66, label = 'DOPC 30b Run 2 (Chain B)')
plt.axvline(O40b1chB[0], 0, 1, linewidth = 2, color = 'yellow', linestyle = '-.', label = 'DOPC 40b Run 1 (Chain B)')
plt.axvline(O40b2chB[0], 0, 1, linewidth = 2, color = 'yellow', linestyle = '-.', alpha = 0.66, label = 'DOPC 40b Run 2 (Chain B)')
plt.axvline(O40b3chB[0], 0, 1, linewidth = 2, color = 'yellow', linestyle = '-.', alpha = 0.33, label = 'DOPC 1b Run 3 (Chain B)')
plt.axvline(O50b1chB[0], 0, 1, linewidth = 2, color = 'yellow', linestyle = ':', label = 'DOPC 50b Run 1 (Chain B)')
plt.axvline(O50b2chB[0], 0, 1, linewidth = 2, color = 'yellow', linestyle = ':', alpha = 0.66, label = 'DOPC 50b Run 2')
plt.axvline(V1b1chB[0], 0, 1, linewidth = 2, color = 'magenta', label = 'DVPC 1b Run 1 (Chain B)')
plt.axvline(V1b2chB[0], 0, 1, linewidth = 2, color = 'magenta', alpha = 0.66, label = 'DVPC 1b Run 2 (Chain B)')
plt.axvline(V40b1chB[0], 0, 1, linewidth = 2, color = 'magenta', linestyle = '--', label = 'DVPC 40b Run 1 (Chain B)')
plt.axvline(V40b2chB[0], 0, 1, linewidth = 2, color = 'magenta', linestyle = '--', alpha = 0.66, label = 'DVPC 40b Run 2 (Chain B)')
plt.axvline(V40b3chB[0], 0, 1, linewidth = 2, color = 'magenta', linestyle = '--', alpha = 0.33, label = 'DVPC 40b Run 3 (Chain B)')
plt.axvline(V50b1chB[0], 0, 1, linewidth = 2, color = 'magenta', linestyle = '-.', label = 'DVPC 50b Run 1 (Chain B)')
plt.axvline(V50b2chB[0], 0, 1, linewidth = 2, color = 'magenta', linestyle = '-.', alpha = 0.66, label = 'DVPC 50b Run 2 (Chain B)')
plt.axvline(V60b1chB[0], 0, 1, linewidth = 2, color = 'magenta', linestyle = ':', label = 'DVPC 60b Run 1 (Chain B)')
plt.axvline(V60b2chB[0], 0, 1, linewidth = 2, color = 'magenta', linestyle = ':', alpha = 0.66, label = 'DVPC 60b Run 2 (Chain B)')
plt.axvline(P1b1chB[0], 0, 1, linewidth = 2, color = 'grey', label = 'POPC 1b Run 1 (Chain B)')
plt.axvline(P1b2chB[0], 0, 1, linewidth = 2, color = 'grey', alpha = 0.66, label = 'POPC 1b Run 2 (Chain B)')
plt.axvline(P50b1chB[0], 0, 1, linewidth = 2, color = 'grey', linestyle = '--', label = 'POPC 50b Run 1 (Chain B)')
plt.axvline(P50b2chB[0], 0, 1, linewidth = 2, color = 'grey', linestyle = '--', alpha = 0.66, label = 'POPC 50b Run 2 (Chain B)')

v = [0,100,0,1]
plt.axis(v)
plt.title("Time when F < 5.5, E > 11.5, Z > 11.5 ($\AA$)", fontsize = 24)
plt.xlabel("Time (ns)", fontsize = '24')
plt.setp(ax.get_yticklabels(), visible = False)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
ax.legend(fontsize = 'small', loc = 'center left', bbox_to_anchor=(1.0,0.5), ncol = 2)
plt.savefig("FEZ-overall-11-11-5.png", format='png', dpi=300)
plt.savefig("FEZ-overall-11-11-5.svg", format='svg', dpi=300)

plt.clf()

dfpc = []
dfpc.append([F1b1chA[0], F1b1chB[0], F1b2chA[0], F1b2chB[0], F30b1chA[0], F30b1chB[0], F30b2chA[0], F30b2chB[0], F40b1chA[0], F40b1chB[0], F40b2chA[0], F40b2chB[0], F40b3chA[0], F40b3chB[0]])

dopc = []
dopc.append([O1b1chA[0], O1b1chB[0], O1b2chA[0], O1b2chB[0], O30b1chA[0], O30b1chB[0], O30b2chA[0], O30b2chB[0], O40b1chA[0], O40b1chB[0], O40b2chA[0], O40b2chB[0], O40b3chA[0], O40b3chB[0], O50b1chA[0], O50b1chB[0], O50b2chA[0], O50b2chB[0]])

dvpc = []
dvpc.append([V1b1chA[0], V1b1chB[0], V1b2chA[0], V1b2chB[0], V40b1chA[0], V40b1chB[0], V40b2chA[0], V40b2chB[0], V40b3chA[0], V40b3chB[0], V50b1chA[0], V50b1chB[0], V50b2chA[0], V50b2chB[0], V60b1chA[0], V60b1chB[0], V60b2chA[0], V60b2chB[0]])

popc = []
popc.append([P1b1chA[0], P1b1chB[0], P1b2chA[0], P1b2chB[0], P50b1chA[0], P50b1chB[0], P50b2chA[0], P50b2chB[0]])

dfpc = np.squeeze(dfpc)
s = pd.Series(dfpc, name = '16:1 $\Delta$9-cis')

dopc = np.squeeze(dopc)
t = pd.Series(dopc, name = '18:1 $\Delta$9-cis')

dvpc = np.squeeze(dvpc)
u = pd.Series(dvpc, name = '20:1 $\Delta$9-cis')

popc = np.squeeze(popc)
w = pd.Series(popc, name = 'POPC')

d = {'16:1 $\Delta$9-cis': s, '18:1 $\Delta$9-cis': t, '20:1 $\Delta$9-cis': u, 'POPC': w}
df = pd.DataFrame(d)

D = {'16:1 $\Delta$9-cis': s, '18:1 $\Delta$9-cis': t, '20:1 $\Delta$9-cis': u}
DF = pd.DataFrame(D)

h = [0,10,0,100]

bx = plt.figure(figsize=(8,8))
sns.boxplot(data=df, whis = 100)
plt.ylabel('Time of state transition (ns)')
plt.xlabel('Lipid Species')
plt.axis(h)
plt.title("Time when F < 5.5, E > 11.5, Z > 11.5 ($\AA$)", fontsize = 20)
plt.savefig("fez-whiskers.png", format='png', dpi=300)
plt.savefig("fez-whiskers.svg", format='svg', dpi=300) 

plt.clf()
plt.figure(figsize=(8,8))
sns.violinplot(data=df, cut = 0)
plt.ylabel('Time of state transition (ns)')
plt.xlabel('Lipid Species')
plt.axis(h)
plt.title("Time when F < 5.5, E > 11.5, Z > 11.5 ($\AA$)", fontsize = 20)
plt.savefig("fez-violin.png", format='png', dpi=300)
plt.savefig("fez-violin.svg", format='svg', dpi=300) 

plt.clf()
sns.boxplot(data=DF, whis = 100)
plt.figure(figsize=(8,8))
plt.ylabel('Time of state transition (ns)')
plt.xlabel('Lipid Species')
plt.axis(h)
plt.title("Time when F < 5.5, E > 11.5, Z > 11.5 ($\AA$)", fontsize = 20)
plt.savefig("fez-whiskers-nopop.png", format='png', dpi=300)
plt.savefig("fez-whiskers-nopop.svg", format='svg', dpi=300) 

plt.clf()
plt.figure(figsize=(8,8))
sns.violinplot(data=DF, cut = 0)
plt.ylabel('Time of state transition (ns)')
plt.xlabel('Lipid Species')
plt.axis(h)
plt.title("Time when F < 5.5, E > 11.5, Z > 11.5 ($\AA$)", fontsize = 20)
plt.savefig("fez-violin-nopop.png", format='png', dpi=300)
plt.savefig("fez-violin-nopop.svg", format='svg', dpi=300) 
