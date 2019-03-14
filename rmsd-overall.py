import numpy as np
import MDAnalysis
import MDAnalysis.analysis.rms
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import pandas as pd
import seaborn as sns

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False


##########
# IMPORT UP STRUCTURE
##########

up = MDAnalysis.Universe('../../../../USEFUL/structures/UP_correctpro.pdb')
prop = up.select_atoms('protein')
selp = prop.select_atoms('(resid 72-98 or resid 183-229 or resid 232-264 or resid 299-330) and name CA')
two = selp.positions.copy()

##########
# DEFINITIONS
##########

def RMSD_traj(A, B, C, D) :
	u = MDAnalysis.Universe(A, B)
	for ts in u.trajectory:
		protein = u.select_atoms('protein')
		sel1 = protein.select_atoms('(resid 72-98 or resid 183-229 or resid 232-264 or resid 299-330) and name CA')
		one = sel1.positions.copy()
		C.append(MDAnalysis.analysis.rms.rmsd(one, two, superposition = True))
	for i in C :
		if i < 2.1 :
			D.append(C.index(i) + 1)
	if D == [] :
		D.append(float(101))

##########
# CALCULATIONS
##########

##### DFPC
### 1b

F1b1RMS = []
F1b1TS = []
RMSD_traj('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc', F1b1RMS, F1b1TS)

F1b2RMS = []
F1b2TS = []
RMSD_traj('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', F1b2RMS, F1b2TS)

### -30b

F30b1RMS = []
F30b1TS = []
RMSD_traj('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', F30b1RMS, F30b1TS)

F30b2RMS = []
F30b2TS = []
RMSD_traj('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', F30b2RMS, F30b2TS)

### -40b

F40b1RMS = []
F40b1TS = []
RMSD_traj('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', F40b1RMS, F40b1TS)

F40b2RMS = []
F40b2TS = []
RMSD_traj('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', F40b2RMS, F40b2TS)

F40b3RMS = []
F40b3TS = []
RMSD_traj('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', F40b3RMS, F40b3TS)

##### DOPC
### 1b

O1b1RMS = []
O1b1TS = []
RMSD_traj('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc', O1b1RMS, O1b1TS)

O1b2RMS = []
O1b2TS = []
RMSD_traj('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', O1b1RMS, O1b2TS)

### -30b

O30b1RMS = []
O30b1TS = []
RMSD_traj('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', O30b1RMS, O30b1TS)

O30b2RMS = []
O30b2TS = []
RMSD_traj('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', O30b2RMS, O30b2TS)

### -40b

O40b1RMS = []
O40b1TS = []
RMSD_traj('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', O40b1RMS, O40b1TS)

O40b2RMS = []
O40b2TS = []
RMSD_traj('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', O40b2RMS, O40b2TS)

O40b3RMS = []
O40b3TS = []
RMSD_traj('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', O40b3RMS, O40b3TS)

### -50b

O50b1RMS = []
O50b1TS = []
RMSD_traj('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', O50b1RMS, O50b1TS)

O50b2RMS = []
O50b2TS = []
RMSD_traj('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', O50b2RMS, O50b2TS)

##### DVPC
### 1b

V1b1RMS = []
V1b1TS = []
RMSD_traj('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc', V1b1RMS, V1b1TS)

V1b2RMS = []
V1b2TS = []
RMSD_traj('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', V1b2RMS, V1b2TS)

### -40b

V40b1RMS = []
V40b1TS = []
RMSD_traj('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', V40b1RMS, V40b1TS)

V40b2RMS = []
V40b2TS = []
RMSD_traj('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', V40b2RMS, V40b2TS)

V40b3RMS = []
V40b3TS = []
RMSD_traj('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', V40b3RMS, V40b3TS)

### -50b

V50b1RMS = []
V50b1TS = []
RMSD_traj('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', V50b1RMS, V50b1TS)

V50b2RMS = []
V50b2TS = []
RMSD_traj('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', V50b2RMS, V50b2TS)

### -60b

V60b1RMS = []
V60b1TS = []
RMSD_traj('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', V60b1RMS, V60b1TS)

V60b2RMS = []
V60b2TS = []
RMSD_traj('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', V60b2RMS, V60b2TS)

### POPC

### 1b

P1b1RMS = []
P1b1TS = []
RMSD_traj('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc', P1b1RMS, P1b1TS)

P1b2RMS = []
P1b2TS = []
RMSD_traj('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc', P1b2RMS, P1b2TS)

### -50b

P50b1RMS = []
P50b1TS = []
RMSD_traj('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', P50b1RMS, P50b1TS)

P50b2RMS = []
P50b2TS = []
RMSD_traj('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', P50b2RMS, P50b2TS)

##########
# PLOTTING
##########

fig = plt.figure(figsize=(16,6))
ax = plt.subplot(111)
plt.axvline(F1b1TS[0], 0, 1, linewidth = 2, color = 'blue', label = 'DFPC 1b Run 1')
plt.axvline(F1b2TS[0], 0, 1, linewidth = 2, color = 'blue', alpha = 0.66, label = 'DFPC 1b Run 2')
plt.axvline(F30b1TS[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '--', label = 'DFPC 30b Run 1')
plt.axvline(F30b2TS[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '--', alpha = 0.66, label = 'DFPC 30b Run 2')
plt.axvline(F40b1TS[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '-.', label = 'DFPC 40b Run 1')
plt.axvline(F40b2TS[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '-.', alpha = 0.66, label = 'DFPC 40b Run 2')
plt.axvline(F40b3TS[0], 0, 1, linewidth = 2, color = 'blue', linestyle = '-.', alpha = 0.33, label = 'DFPC 40b Run 3')
plt.axvline(O1b1TS[0], 0, 1, linewidth = 2, color = 'green', label = 'DOPC 1b Run 1')
plt.axvline(O1b2TS[0], 0, 1, linewidth = 2, color = 'green', alpha = 0.66, label = 'DOPC 1b Run 2')
plt.axvline(O30b1TS[0], 0, 1, linewidth = 2, color = 'green', linestyle = '--', label = 'DOPC 30b Run 1')
plt.axvline(O30b2TS[0], 0, 1, linewidth = 2, color = 'green', linestyle = '--', alpha = 0.66, label = 'DOPC 30b Run 2')
plt.axvline(O40b1TS[0], 0, 1, linewidth = 2, color = 'green', linestyle = '-.', label = 'DOPC 40b Run 1')
plt.axvline(O40b2TS[0], 0, 1, linewidth = 2, color = 'green', linestyle = '-.', alpha = 0.66, label = 'DOPC 40b Run 2')
plt.axvline(O40b3TS[0], 0, 1, linewidth = 2, color = 'green', linestyle = '-.', alpha = 0.33, label = 'DOPC 1b Run 3')
plt.axvline(O50b1TS[0], 0, 1, linewidth = 2, color = 'green', linestyle = ':', label = 'DOPC 50b Run 1')
plt.axvline(O50b2TS[0], 0, 1, linewidth = 2, color = 'green', linestyle = ':', alpha = 0.66, label = 'DOPC 50b Run 2')
plt.axvline(V1b1TS[0], 0, 1, linewidth = 2, color = 'red', label = 'DVPC 1b Run 1')
plt.axvline(V1b2TS[0], 0, 1, linewidth = 2, color = 'red', alpha = 0.66, label = 'DVPC 1b Run 2')
plt.axvline(V40b1TS[0], 0, 1, linewidth = 2, color = 'red', linestyle = '--', label = 'DVPC 40b Run 1')
plt.axvline(V40b2TS[0], 0, 1, linewidth = 2, color = 'red', linestyle = '--', alpha = 0.66, label = 'DVPC 40b Run 2')
plt.axvline(V40b3TS[0], 0, 1, linewidth = 2, color = 'red', linestyle = '--', alpha = 0.33, label = 'DVPC 40b Run 3')
plt.axvline(V50b1TS[0], 0, 1, linewidth = 2, color = 'red', linestyle = '-.', label = 'DVPC 50b Run 1')
plt.axvline(V50b2TS[0], 0, 1, linewidth = 2, color = 'red', linestyle = '-.', alpha = 0.66, label = 'DVPC 50b Run 2')
plt.axvline(V60b1TS[0], 0, 1, linewidth = 2, color = 'red', linestyle = ':', label = 'DVPC 60b Run 1')
plt.axvline(V60b2TS[0], 0, 1, linewidth = 2, color = 'red', linestyle = ':', alpha = 0.66, label = 'DVPC 60b Run 2')
plt.axvline(P1b1TS[0], 0, 1, linewidth = 2, color = 'black', linestyle = '-.', label = 'POPC 1b Run 1')
plt.axvline(P1b2TS[0], 0, 1, linewidth = 2, color = 'black', linestyle = '-.', alpha = 0.66, label = 'POPC 1b Run 2')
plt.axvline(P50b1TS[0], 0, 1, linewidth = 2, color = 'black', linestyle = ':', label = 'POPC 50b Run 1')
plt.axvline(P50b2TS[0], 0, 1, linewidth = 2, color = 'black', linestyle = ':', alpha = 0.66, label = 'POPC 50b Run 2')
v = [0,100,0,1]
plt.axis(v)
plt.title("Time when RMSD < 2.1 ($\AA$)", fontsize = 24)
plt.xlabel("Time (ns)", fontsize = '24')
plt.setp(ax.get_yticklabels(), visible = False)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(fontsize = 'medium', loc = 'center left', bbox_to_anchor=(1.05,0.5), ncol = 2)
plt.savefig("rmsd-overall.png", format='png', dpi=300)
plt.savefig("rmsd-overall.svg", format='svg', dpi=300)

plt.clf()

dfpc = []
dfpc.append([F1b1TS[0], F1b2TS[0], F30b1TS[0], F30b2TS[0], F40b1TS[0], F40b2TS[0], F40b3TS[0]])

dopc = []
dopc.append([O1b1TS[0], O1b2TS[0], O30b1TS[0], O30b2TS[0], O40b1TS[0], O40b2TS[0], O40b3TS[0], O50b1TS[0], O50b2TS[0]])

dvpc = []
dvpc.append([V1b1TS[0], V1b2TS[0], V40b1TS[0], V40b2TS[0], V40b3TS[0], V50b1TS[0], V50b2TS[0], V60b1TS[0], V60b2TS[0]])

popc = []
popc.append([P1b1TS[0], P1b2TS[0], P50b1TS[0], P50b2TS[0]])

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

bx = plt.figure(figsize=(8,8))
sns.boxplot(data=df, whis = 100)
plt.ylabel('Time (ns)')
plt.xlabel('Lipid Species')
plt.title("Time when RMSD < 2.1 ($\AA$)", fontsize = 20)
plt.savefig("rmsd-whiskers.png", format='png', dpi=300)
plt.savefig("rmsd-whiskers.svg", format='svg', dpi=300) 

plt.clf()
plt.figure(figsize=(8,8))
sns.violinplot(data=df, cut = 0)
plt.ylabel('Time (ns)')
plt.xlabel('Lipid Species')
plt.title("Time when RMSD < 2.1 ($\AA$)", fontsize = 20)
plt.savefig("rmsd-violin.png", format='png', dpi=300)
plt.savefig("rmsd-violin.svg", format='svg', dpi=300) 

bx = plt.figure(figsize=(8,8))
sns.boxplot(data=DF, whis = 100)
plt.ylabel('Time (ns)')
plt.xlabel('Lipid Species')
plt.title("Time when RMSD < 2.1 ($\AA$)", fontsize = 20)
plt.savefig("rmsd-whiskers-nopop.png", format='png', dpi=300)
plt.savefig("rmsd-whiskers-nopop.svg", format='svg', dpi=300) 

plt.clf()
plt.figure(figsize=(8,8))
sns.violinplot(data=DF, cut = 0)
plt.ylabel('Time (ns)')
plt.xlabel('Lipid Species')
plt.title("Time when RMSD < 2.1 ($\AA$)", fontsize = 20)
plt.savefig("rmsd-violin-nopop.png", format='png', dpi=300)
plt.savefig("rmsd-violin-nopop.svg", format='svg', dpi=300) 

