import MDAnalysis
import numpy as np
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
# DEFINITIONS
##########

def dif (A):
    u = MDAnalysis.Universe(A)
    protein = u.select_atoms('protein')
    L = u.select_atoms('name P*')
    proCOM = protein.center_of_mass()
    memCOM = L.center_of_mass()
    dif = proCOM[2] - memCOM[2]
    return dif

def meandif (A, B):
    u = MDAnalysis.Universe(A, B)
    protein = u.select_atoms('protein')
    L = u.select_atoms('name P*')
    meanrun = []
    for ts in u.trajectory:
        proCOM = protein.center_of_mass()
        memCOM = L.center_of_mass()
        meanrun.append(proCOM[2] - memCOM[2])
    dif = np.mean(meanrun[25:75])
    return dif

##########
# CALCULATIONS
##########

Dif = []
Meandif = []

##### DFPC
### 1b

Dif.append(dif('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb'))
Meandif.append(meandif('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc'))

Dif.append(dif('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb'))
Meandif.append(meandif('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc'))

### -30b

Dif.append(dif('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc'))

### -40b

Dif.append(dif('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc'))

Dif.append(dif('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb'))
Meandif.append(meandif('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc'))

##### DOPC
### 1b

Dif.append(dif('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc'))

Dif.append(dif('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc'))

### -30b

Dif.append(dif('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc'))

### -40b

Dif.append(dif('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc'))

Dif.append(dif('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc'))

### -50b

Dif.append(dif('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc'))

##### DVPC
### 1b

Dif.append(dif('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc'))

Dif.append(dif('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc'))

### -40b

Dif.append(dif('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc'))

Dif.append(dif('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc'))

### -50b

Dif.append(dif('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc'))

### -60b

Dif.append(dif('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc'))

Dif.append(dif('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc'))

##### POPC
### 1b

Dif.append(dif('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc'))

### -50b

Dif.append(dif('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb'))
Meandif.append(meandif('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc'))

Dif.append(dif('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb'))
Meandif.append(meandif('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc'))

###

uni = MDAnalysis.Universe('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb')
prot = uni.select_atoms('protein')
Lipid = uni.select_atoms('name P*')
membCOM = Lipid.center_of_mass()
protCOM = prot.center_of_mass()
calphas = prot.select_atoms('name CA').positions

##########
# PLOTTING
##########

plt.figure(figsize=(14,10))

plt.subplot(121)
plt.plot(calphas[:,0], calphas[:,2], color = 'purple', linewidth = 2)
plt.axhline(protCOM[2], color = 'black', linestyle = '--')
plt.axhline(Dif[0] + protCOM[2], label = 'DFPC Start', color = 'blue')
plt.axhline(Dif[7] + protCOM[2], label = 'DOPC Start', color = 'green')
plt.axhline(Dif[16] + protCOM[2], label = 'DVPC Start', color = 'red')
plt.axhline(Dif[25] + protCOM[2], label = 'POPC Start', color = 'black')

plt.axhline(Meandif[0] + protCOM[2], label = 'F1b1 Average', color = 'blue', linewidth = 1)
plt.axhline(Meandif[1] + protCOM[2], label = 'F1b2 Average', color = 'blue', linewidth = 1, alpha = 0.66)
plt.axhline(Meandif[2] + protCOM[2], label = 'F30b1 Average', color = 'blue', linewidth = 1, linestyle = '--')
plt.axhline(Meandif[3] + protCOM[2], label = 'F30b2 Average', color = 'blue', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline(Meandif[4] + protCOM[2], label = 'F40b1 Average', color = 'blue', linewidth = 1, linestyle = '-.')
plt.axhline(Meandif[5] + protCOM[2], label = 'F40b2 Average', color = 'blue', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline(Meandif[6] + protCOM[2], label = 'F40b3 Average', color = 'blue', linewidth = 1, linestyle = '-.', alpha = 0.33)

plt.axhline(Meandif[7] + protCOM[2], label = 'O1b1 Average', color = 'green', linewidth = 1)
plt.axhline(Meandif[8] + protCOM[2], label = 'O1b2 Average', color = 'green', linewidth = 1, alpha = 0.66)
plt.axhline(Meandif[9] + protCOM[2], label = 'O30b1 Average', color = 'green', linewidth = 1, linestyle = '--')
plt.axhline(Meandif[10] + protCOM[2], label = 'O30b2 Average', color = 'green', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline(Meandif[11] + protCOM[2], label = 'O40b1 Average', color = 'green', linewidth = 1, linestyle = '-.')
plt.axhline(Meandif[12] + protCOM[2], label = 'O40b2 Average', color = 'green', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline(Meandif[13] + protCOM[2], label = 'O40b3 Average', color = 'green', linewidth = 1, linestyle = '-.', alpha = 0.33)
plt.axhline(Meandif[14] + protCOM[2], label = 'O50b1 Average', color = 'green', linewidth = 1, linestyle = ':')
plt.axhline(Meandif[15] + protCOM[2], label = 'O50b2 Average', color = 'green', linewidth = 1, linestyle = ':', alpha = 0.66)

plt.axhline(Meandif[16] + protCOM[2], label = 'V1b1 Average', color = 'red', linewidth = 1)
plt.axhline(Meandif[17] + protCOM[2], label = 'V1b2 Average', color = 'red', linewidth = 1, alpha = 0.66)
plt.axhline(Meandif[18] + protCOM[2], label = 'V40b1 Average', color = 'red', linewidth = 1, linestyle = '--')
plt.axhline(Meandif[19] + protCOM[2], label = 'V40b2 Average', color = 'red', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline(Meandif[20] + protCOM[2], label = 'V40b3 Average', color = 'red', linewidth = 1, linestyle = '--', alpha = 0.33)
plt.axhline(Meandif[21] + protCOM[2], label = 'V50b1 Average', color = 'red', linewidth = 1, linestyle = '-.')
plt.axhline(Meandif[22] + protCOM[2], label = 'V50b2 Average', color = 'red', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline(Meandif[23] + protCOM[2], label = 'V60b1 Average', color = 'red', linewidth = 1, linestyle = ':')
plt.axhline(Meandif[24] + protCOM[2], label = 'V60b2 Average', color = 'red', linewidth = 1, linestyle = ':', alpha = 0.66)

plt.axhline(Meandif[25] + protCOM[2], label = 'P1b1 Average', color = 'black', linewidth = 1)
plt.axhline(Meandif[26] + protCOM[2], label = 'P1b2 Average', color = 'black', linewidth = 1, alpha = 0.66)
plt.axhline(Meandif[27] + protCOM[2], label = 'P50b1 Average', color = 'black', linewidth = 1, linestyle = '--')
plt.axhline(Meandif[28] + protCOM[2], label = 'P50b2 Average', color = 'black', linewidth = 1, linestyle = '--', alpha = 0.66)

plt.ylabel("z coordinate ($\AA$)", fontsize = '18')
plt.xlabel("x coordinate ($\AA$)", fontsize = '18')
w = [20,90,20,120]
plt.axis(w)

ax = plt.subplot(122)

plt.axhline(Dif[0], label = 'DFPC Start', color = 'blue', linewidth = 3)
plt.axhline(Dif[7], label = 'DOPC Start', color = 'green', linewidth = 3)
plt.axhline(Dif[16], label = 'DVPC Start', color = 'red', linewidth = 3)
plt.axhline(Dif[25], label = 'POPC Start', color = 'black', linewidth = 3)

plt.axhline(Meandif[0], label = 'F1b1 Average', color = 'blue', linewidth = 2)
plt.axhline(Meandif[1], label = 'F1b2 Average', color = 'blue', linewidth = 2, alpha = 0.66)
plt.axhline(Meandif[2], label = 'F30b1 Average', color = 'blue', linewidth = 2, linestyle = '--')
plt.axhline(Meandif[3], label = 'F30b2 Average', color = 'blue', linewidth = 2, linestyle = '--', alpha = 0.66)
plt.axhline(Meandif[4], label = 'F40b1 Average', color = 'blue', linewidth = 2, linestyle = '-.')
plt.axhline(Meandif[5], label = 'F40b2 Average', color = 'blue', linewidth = 2, linestyle = '-.', alpha = 0.66)
plt.axhline(Meandif[6], label = 'F40b3 Average', color = 'blue', linewidth = 2, linestyle = '-.', alpha = 0.33)

plt.axhline(Meandif[7], label = 'O1b1 Average', color = 'green', linewidth = 2)
plt.axhline(Meandif[8], label = 'O1b2 Average', color = 'green', linewidth = 2, alpha = 0.66)
plt.axhline(Meandif[9], label = 'O30b1 Average', color = 'green', linewidth = 2, linestyle = '--')
plt.axhline(Meandif[10], label = 'O30b2 Average', color = 'green', linewidth = 2, linestyle = '--', alpha = 0.66)
plt.axhline(Meandif[11], label = 'O40b1 Average', color = 'green', linewidth = 2, linestyle = '-.')
plt.axhline(Meandif[12], label = 'O40b2 Average', color = 'green', linewidth = 2, linestyle = '-.', alpha = 0.66)
plt.axhline(Meandif[13], label = 'O40b3 Average', color = 'green', linewidth = 2, linestyle = '-.', alpha = 0.33)
plt.axhline(Meandif[14], label = 'O50b1 Average', color = 'green', linewidth = 2, linestyle = ':')
plt.axhline(Meandif[15], label = 'O50b2 Average', color = 'green', linewidth = 2, linestyle = ':', alpha = 0.66)

plt.axhline(Meandif[16], label = 'V1b1 Average', color = 'red', linewidth = 2)
plt.axhline(Meandif[17], label = 'V1b2 Average', color = 'red', linewidth = 2, alpha = 0.66)
plt.axhline(Meandif[18], label = 'V40b1 Average', color = 'red', linewidth = 2, linestyle = '--')
plt.axhline(Meandif[19], label = 'V40b2 Average', color = 'red', linewidth = 2, linestyle = '--', alpha = 0.66)
plt.axhline(Meandif[20], label = 'V40b3 Average', color = 'red', linewidth = 2, linestyle = '--', alpha = 0.33)
plt.axhline(Meandif[21], label = 'V50b1 Average', color = 'red', linewidth = 2, linestyle = '-.')
plt.axhline(Meandif[22], label = 'V50b2 Average', color = 'red', linewidth = 2, linestyle = '-.', alpha = 0.66)
plt.axhline(Meandif[23], label = 'V60b1 Average', color = 'red', linewidth = 2, linestyle = ':')
plt.axhline(Meandif[24], label = 'V60b2 Average', color = 'red', linewidth = 2, linestyle = ':', alpha = 0.66)

plt.axhline(Meandif[25], label = 'P1b1 Average', color = 'black', linewidth = 2)
plt.axhline(Meandif[26], label = 'P1b2 Average', color = 'black', linewidth = 2, alpha = 0.66)
plt.axhline(Meandif[27], label = 'P50b1 Average', color = 'black', linewidth = 2, linestyle = '--')
plt.axhline(Meandif[28], label = 'P50b2 Average', color = 'black', linewidth = 2, linestyle = '--', alpha = 0.66)

plt.ylabel("Distance ($\AA$)", fontsize = '18')
plt.title('Height of Protein Centre-of-Mass\nabove Bilayer Centre-of-Mass (0 $\AA$)', fontsize = 24)
v = [0,1,55.5,58.5]
plt.axis(v)
plt.setp(ax.get_xticklabels(), visible = False)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(fontsize = 'large', loc = 'center left', bbox_to_anchor=(1.05,0.5))

plt.suptitle('Position of bilayer centre-of-mass', fontsize = 36)
plt.savefig("proposoverall.png", format='png', dpi=300)
plt.savefig("proposoverall.svg", format='svg', dpi=300)

plt.clf()

dfpc = []
dfpc.append([Meandif[0], Meandif[1], Meandif[2], Meandif[3], Meandif[4], Meandif[5], Meandif[6]])

dopc = []
dopc.append([Meandif[7], Meandif[8], Meandif[9], Meandif[10], Meandif[11], Meandif[12], Meandif[13], Meandif[14], Meandif[15]])

dvpc = []
dvpc.append([Meandif[16], Meandif[17], Meandif[18], Meandif[19], Meandif[20], Meandif[21], Meandif[22], Meandif[23], Meandif[24]])

popc = []
popc.append([Meandif[25], Meandif[26], Meandif[27], Meandif[28]])

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

#bx = plt.figure(figsize=(8,8))
#sns.boxplot(data=df, whis = 100)
#plt.ylabel('Distance ($\AA$)')
#plt.xlabel('Lipid Species')
#plt.title('Height of Protein Centre-of-Mass\nabove Bilayer Centre-of-Mass (0 $\AA$)', fontsize = 24)
#plt.savefig("propos-whiskers.png", format='png', dpi=300)
#plt.savefig("propos-whiskers.svg", format='svg', dpi=300) 

#plt.clf()
#plt.figure(figsize=(8,8))
#sns.violinplot(data=df, cut = 0)
#plt.ylabel('Distance ($\AA$)')
#plt.xlabel('Lipid Species')
#plt.title('Height of Protein Centre-of-Mass\nabove Bilayer Centre-of-Mass (0 $\AA$)', fontsize = 24)
#plt.savefig("propos-violin.png", format='png', dpi=300)
#plt.savefig("propos-violin.svg", format='svg', dpi=300) 

#plt.clf()
#plt.figure(figsize=(8,8))
#sns.boxplot(data=DF, whis = 100)
#plt.ylabel('Distance ($\AA$)')
#plt.xlabel('Lipid Species')
#plt.title('Height of Protein Centre-of-Mass\nabove Bilayer Centre-of-Mass (0 $\AA$)', fontsize = 24)
#plt.savefig("propos-whiskers-nopop.png", format='png', dpi=300)
#plt.savefig("propos-whiskers-nopop.svg", format='svg', dpi=300) 

#plt.clf()
#plt.figure(figsize=(8,8))
#sns.violinplot(data=DF, cut = 0)
#plt.ylabel('Distance ($\AA$)')
#plt.xlabel('Lipid Species')
#plt.title('Height of Protein Centre-of-Mass\nabove Bilayer Centre-of-Mass (0 $\AA$)', fontsize = 24)
#plt.savefig("propos-violin-nopop.png", format='png', dpi=300)
#plt.savefig("propos-violin-nopop.svg", format='svg', dpi=300) 

plt.figure(figsize=(6,8))
sns.stripplot(data=DF, jitter = 0, size = 10, linewidth = 1)
plt.ylabel('Distance ($\AA$)')
plt.xlabel('Lipid Species')
#plt.axis(h)
plt.title("Height of Protein Centre-of-Mass\nabove Bilayer Centre-of-Mass (0 $\AA$)", fontsize = 24)
plt.savefig("propos-strip-nopop.png", format='png', dpi=300)
plt.savefig("propos-strip-nopop.svg", format='svg', dpi=300) 
