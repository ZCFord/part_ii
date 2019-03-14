import numpy as np
import MDAnalysis
import MDAnalysis.analysis.leaflet
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import seaborn as sns
import pandas as pd

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False

##########
# CALCULATIONS
##########

##### DFPC
### 1b

d1 = np.loadtxt('../../DFPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)
mean1 = np.mean(d1[25:75, 1])

d2 = np.loadtxt('../../DFPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)
mean2 = np.mean(d2[25:75, 1])

av1 = (mean1+mean2)/2

### -30b

d3 = np.loadtxt('../../DFPC/Backwards/md/30bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean3 = np.mean(d3[25:75, 1])

d4 = np.loadtxt('../../DFPC/Backwards/md/30bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean4 = np.mean(d4[25:75, 1])

av2 = (mean3+mean4)/2

### -40b

d5 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean5 = np.mean(d5[25:75, 1])

d6 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean6 = np.mean(d6[25:75, 1])

d7 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)
mean7 = np.mean(d7[25:75, 1])

av3 = (mean5+mean6+mean7)/3

##### DOPC
### 1b

d8 = np.loadtxt('../../DOPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)
mean8 = np.mean(d8[25:75, 1])

d9 = np.loadtxt('../../DOPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)
mean9 = np.mean(d9[25:75, 1])

av4 = (mean8+mean9)/2

### -30b

d10 = np.loadtxt('../../DOPC/Backwards/md/30bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean10 = np.mean(d10[25:75, 1])

d11 = np.loadtxt('../../DOPC/Backwards/md/30bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean11 = np.mean(d11[25:75, 1])

av5 = (mean10+mean11)/2

### -40b

d12 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean12 = np.mean(d12[25:75, 1])

d13 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean13 = np.mean(d13[25:75, 1])

d14 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)
mean14 = np.mean(d14[25:75, 1])

av6 = (mean12+mean13+mean14)/3

### -50b

d15 = np.loadtxt('../../DOPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean15 = np.mean(d15[25:75, 1])

d16 = np.loadtxt('../../DOPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean16 = np.mean(d16[25:75, 1])

av7 = (mean15+mean16)/2

##### DVPC
### 1b

d17 = np.loadtxt('../../DVPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)
mean17 = np.mean(d17[25:75, 1])

d18 = np.loadtxt('../../DVPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)
mean18 = np.mean(d18[25:75, 1])

av8 = (mean17+mean18)/2

### -40b

d19 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean19 = np.mean(d19[25:75, 1])

d20 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean20 = np.mean(d20[25:75, 1])

d21 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)
mean21 = np.mean(d21[25:75, 1])

av9 = (mean19+mean20+mean21)/3

### -50b

d22 = np.loadtxt('../../DVPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean22 = np.mean(d22[25:75, 1])

d23 = np.loadtxt('../../DVPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean23 = np.mean(d23[25:75, 1])

av10 = (mean22+mean23)/2

### -60b

d24 = np.loadtxt('../../DVPC/Backwards/md/60bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean24 = np.mean(d24[25:75, 1])

d25 = np.loadtxt('../../DVPC/Backwards/md/60bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean25 = np.mean(d25[25:75, 1])

av11 = (mean24+mean25)/2

##### POPC
### 1b

d26 = np.loadtxt('../../POPC/Backwards/md/1bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean26 = np.mean(d26[25:75, 1])

d27 = np.loadtxt('../../POPC/Backwards/md/1bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean27 = np.mean(d27[25:75, 1])

av12 = (mean26+mean27)/2

### -50b

d28 = np.loadtxt('../../POPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean28 = np.mean(d28[25:75, 1])

d29 = np.loadtxt('../../POPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean29 = np.mean(d29[25:75, 1])

av13 = (mean28+mean29)/2

##########
# PLOTTING
##########

fig = plt.figure(figsize=(12,12))
ax = plt.subplot(111)
plt.axhline(av1*100, label = 'DFPC 1b', color = 'blue', linewidth = 2)
plt.axhline(av2*100, label = 'DFPC -30b', color = 'blue', linewidth = 2, linestyle = '--')
plt.axhline(av3*100, label = 'DFPC -40b', color = 'blue', linewidth = 2, linestyle = '-.')
plt.axhline(av4*100, label = 'DOPC 1b', color = 'green', linewidth = 2)
plt.axhline(av5*100, label = 'DOPC -30b', color = 'green', linewidth = 2, linestyle = '--')
plt.axhline(av6*100, label = 'DOPC -40b', color = 'green', linewidth = 2, linestyle = '-.')
plt.axhline(av7*100, label = 'DOPC -50b', color = 'green', linewidth = 2, linestyle = ':')
plt.axhline(av8*100, label = 'DVPC 1b', color = 'red', linewidth = 2)
plt.axhline(av9*100, label = 'DVPC -40b', color = 'red', linewidth = 2, linestyle = '--')
plt.axhline(av10*100, label = 'DVPC -50b', color = 'red', linewidth = 2, linestyle = '-.')
plt.axhline(av11*100, label = 'DVPC -60b', color = 'red', linewidth = 2, linestyle = ':')
plt.axhline(av12*100, label = 'POPC 1b', color = 'black', linewidth = 2)
plt.axhline(av13*100, label = 'POPC -50b', color = 'black', linewidth = 2, linestyle = '--')
v = [0,1,50,120]
plt.axis(v)
plt.ylabel("Area Per Lipid ($\AA^2$)", fontsize = '24')
plt.setp(ax.get_xticklabels(), visible = False)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(fontsize = 'x-large', loc = 'center left', bbox_to_anchor=(1,0.5))
plt.savefig("aploverall.png", format='png', dpi=300)
plt.savefig("aploverall.svg", format='svg', dpi=300)

plt.clf()

dfpc = []
dfpc.append([mean1*100, mean2*100, mean3*100, mean4*100, mean5*100, mean6*100, mean7*100])

dopc = []
dopc.append([mean8*100, mean9*100, mean10*100, mean11*100, mean12*100, mean13*100, mean14*100, mean15*100, mean16*100])

dvpc = []
dvpc.append([mean17*100, mean18*100, mean19*100, mean20*100, mean21*100, mean22*100, mean23*100, mean24*100, mean25*100])

popc = []
popc.append([mean26*100, mean27*100, mean28*100, mean29*100])

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
#plt.ylabel('Area Per Lipid ($\AA^2$)')
#plt.xlabel('Lipid Species')
#plt.savefig("apl-whiskers.png", format='png', dpi=300)
#plt.savefig("apl-whiskers.svg", format='svg', dpi=300) 

#plt.clf()
#plt.figure(figsize=(8,8))
#sns.violinplot(data=df, cut = 0)
#plt.ylabel('Area Per Lipid ($\AA^2$)')
#plt.xlabel('Lipid Species')
#plt.savefig("apl-violin.png", format='png', dpi=300)
#plt.savefig("apl-violin.svg", format='svg', dpi=300) 
    
#plt.clf()          
#plt.figure(figsize=(8,8))
#sns.boxplot(data=DF, whis = 100)
#plt.ylabel('Area Per Lipid ($\AA^2$)')
#plt.xlabel('Lipid Species')
#plt.savefig("apl-whiskers-nopop.png", format='png', dpi=300)
#plt.savefig("apl-whiskers-nopop.svg", format='svg', dpi=300) 

#plt.clf()
#plt.figure(figsize=(8,8))
#sns.violinplot(data=DF, cut = 0)
#plt.ylabel('Area Per Lipid ($\AA^2$)')
#plt.xlabel('Lipid Species')
#plt.savefig("apl-violin-nopop.png", format='png', dpi=300)
#plt.savefig("apl-violin-nopop.svg", format='svg', dpi=300) 
 
plt.figure(figsize=(6,8))
sns.stripplot(data=DF, jitter = 0, size = 10, linewidth = 1)
plt.ylabel('Area Per Lipid ($\AA^2$)')
plt.xlabel('Lipid Species')
plt.savefig("apl-strip-nopop.png", format='png', dpi=300)
plt.savefig("apl-strip-nopop.svg", format='svg', dpi=300)             
