import MDAnalysis
import MDAnalysis.analysis.leaflet
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import pandas as pd
import seaborn as sns
import numpy as np

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False

##########
# DEFINITIONS
##########

def bilayer_thickness (A, B) : # finds the bilayer thickness from a file or universe input
	downP = B.groups(0)
	upP = B.groups(1)
	upPxyz = upP.positions
	upPz = upPxyz[:,2] #selects z coordinates
	upPz_av = np.mean(upPz)
	downPxyz = downP.positions
	downPz = downPxyz[:,2] 
	downPz_av = np.mean(downPz)
	bilayer_thickness = downPz_av - upPz_av #difference in average upper and lower z coordinates
	return (bilayer_thickness)

def dif (A):
    u = MDAnalysis.Universe(A)
    protein = u.select_atoms('protein')
    L = u.select_atoms('name P*')
    proCOM = protein.center_of_mass()
    memCOM = L.center_of_mass()
    dif = memCOM[2] - proCOM[2]
    return dif

def meandif (A, B):
    u = MDAnalysis.Universe(A, B)
    protein = u.select_atoms('protein')
    L = u.select_atoms('name P*')
    meanrun = []
    for ts in u.trajectory:
        proCOM = protein.center_of_mass()
        memCOM = L.center_of_mass()
        meanrun.append(memCOM[2] - proCOM[2])
    dif = np.mean(meanrun[25:75])
    return dif

##########
# CALCULATIONS
##########

# BILAYER THICKNESS

##### DFPC
### 1b

u1 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc')
L1 = MDAnalysis.analysis.leaflet.LeafletFinder(u1, 'name P*')

bilayer1 = []

for ts in u1.trajectory :
	bilayer1.append(bilayer_thickness(u1, L1))

mean1 = np.mean(bilayer1[25:75])

u2 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L2 = MDAnalysis.analysis.leaflet.LeafletFinder(u2, 'name P*')

bilayer2 = []

for ts in u2.trajectory :
	bilayer2.append(bilayer_thickness(u2, L2))

mean2 = np.mean(bilayer2[25:75])

av1 = (mean1+mean2)/2

### -30b

u3 = MDAnalysis.Universe('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc')
L3 = MDAnalysis.analysis.leaflet.LeafletFinder(u3, 'name P*')

bilayer3 = []

for ts in u3.trajectory :
	bilayer3.append(bilayer_thickness(u3, L3))

mean3 = np.mean(bilayer3[25:75])

u4 = MDAnalysis.Universe('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc')
L4 = MDAnalysis.analysis.leaflet.LeafletFinder(u4, 'name P*')

bilayer4 = []

for ts in u4.trajectory :
	bilayer4.append(bilayer_thickness(u4, L4))

mean4 = np.mean(bilayer4[25:75])

av2 = (mean3+mean4)/2

### -40b

bilayer5 = []

u5 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L5 = MDAnalysis.analysis.leaflet.LeafletFinder(u5, 'name P*') #selects P atoms and puts in two leaflet

for ts in u5.trajectory :
	bilayer5.append(bilayer_thickness(u5, L5))

mean5 = np.mean(bilayer5[25:75])

bilayer6 = []

u6 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L6 = MDAnalysis.analysis.leaflet.LeafletFinder(u6, 'name P*') #selects P atoms and puts in two leaflet

for ts in u6.trajectory :
	bilayer6.append(bilayer_thickness(u6, L6))

mean6 = np.mean(bilayer6[25:75])

bilayer7 = []

u7 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L7 = MDAnalysis.analysis.leaflet.LeafletFinder(u7, 'name P*') #selects P atoms and puts in two leaflet

for ts in u7.trajectory :
	bilayer7.append(bilayer_thickness(u7, L7))

mean7 = np.mean(bilayer7[25:75])

av3 = (mean5 + mean6 + mean7)/3

##### DOPC
### 1b

u8 = MDAnalysis.Universe('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc')
L8 = MDAnalysis.analysis.leaflet.LeafletFinder(u8, 'name P*')

bilayer8 = []

for ts in u8.trajectory :
	bilayer8.append(bilayer_thickness(u8, L8))

mean8 = np.mean(bilayer8[25:75])

u9 = MDAnalysis.Universe('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L9 = MDAnalysis.analysis.leaflet.LeafletFinder(u9, 'name P*')

bilayer9 = []

for ts in u9.trajectory :
	bilayer9.append(bilayer_thickness(u9, L9))

mean9 = np.mean(bilayer9[25:75])

av4 = (mean8+mean9)/2

### -30b

u10 = MDAnalysis.Universe('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc')
L10 = MDAnalysis.analysis.leaflet.LeafletFinder(u10, 'name P*')

bilayer10 = []

for ts in u10.trajectory :
	bilayer10.append(bilayer_thickness(u10, L10))

mean10 = np.mean(bilayer10[25:75])

u11 = MDAnalysis.Universe('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc')
L11 = MDAnalysis.analysis.leaflet.LeafletFinder(u11, 'name P*')

bilayer11 = []

for ts in u11.trajectory :
	bilayer11.append(bilayer_thickness(u11, L11))

mean11 = np.mean(bilayer11[25:75])

av5 = (mean10+mean11)/2

### -40b

bilayer12 = []

u12 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L12 = MDAnalysis.analysis.leaflet.LeafletFinder(u12, 'name P*') #selects P atoms and puts in two leaflet

for ts in u12.trajectory :
	bilayer12.append(bilayer_thickness(u12, L12))

mean12 = np.mean(bilayer12[25:75])

bilayer13 = []

u13 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L13 = MDAnalysis.analysis.leaflet.LeafletFinder(u13, 'name P*') #selects P atoms and puts in two leaflet

for ts in u13.trajectory :
	bilayer13.append(bilayer_thickness(u13, L13))

mean13 = np.mean(bilayer13[25:75])

bilayer14 = []

u14 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L14 = MDAnalysis.analysis.leaflet.LeafletFinder(u14, 'name P*') #selects P atoms and puts in two leaflet

for ts in u14.trajectory :
	bilayer14.append(bilayer_thickness(u14, L14))

mean14 = np.mean(bilayer14[25:75])

av6 = (mean12 + mean13 + mean14)/3

### -50b

u15 = MDAnalysis.Universe('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L15 = MDAnalysis.analysis.leaflet.LeafletFinder(u15, 'name P*')

bilayer15 = []

for ts in u15.trajectory :
	bilayer15.append(bilayer_thickness(u15, L15))

mean15 = np.mean(bilayer15[25:75])

u16 = MDAnalysis.Universe('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L16 = MDAnalysis.analysis.leaflet.LeafletFinder(u16, 'name P*')

bilayer16 = []

for ts in u16.trajectory :
	bilayer16.append(bilayer_thickness(u16, L16))

mean16 = np.mean(bilayer16[25:75])

av7 = (mean15+mean16)/2

##### DVPC
### 1b

u17 = MDAnalysis.Universe('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc')
L17 = MDAnalysis.analysis.leaflet.LeafletFinder(u17, 'name P*')

bilayer17 = []

for ts in u17.trajectory :
	bilayer17.append(bilayer_thickness(u17, L17))

mean17 = np.mean(bilayer17[25:75])

u18 = MDAnalysis.Universe('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L18 = MDAnalysis.analysis.leaflet.LeafletFinder(u18, 'name P*')

bilayer18 = []

for ts in u18.trajectory :
	bilayer18.append(bilayer_thickness(u18, L18))

mean18 = np.mean(bilayer18[25:75])

av8 = (mean17+mean18)/2

### -40b

bilayer19 = []

u19 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L19 = MDAnalysis.analysis.leaflet.LeafletFinder(u19, 'name P*') #selects P atoms and puts in two leaflet

for ts in u19.trajectory :
	bilayer19.append(bilayer_thickness(u19, L19))

mean19 = np.mean(bilayer19[25:75])

bilayer20 = []

u20 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L20 = MDAnalysis.analysis.leaflet.LeafletFinder(u20, 'name P*') #selects P atoms and puts in two leaflet

for ts in u20.trajectory :
	bilayer20.append(bilayer_thickness(u20, L20))

mean20 = np.mean(bilayer20[25:75])

bilayer21 = []

u21 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L21 = MDAnalysis.analysis.leaflet.LeafletFinder(u21, 'name P*') #selects P atoms and puts in two leaflet

for ts in u21.trajectory :
	bilayer21.append(bilayer_thickness(u21, L21))

mean21 = np.mean(bilayer21[25:75])

av9 = (mean19 + mean20 + mean21)/3

### -50b

u22 = MDAnalysis.Universe('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L22 = MDAnalysis.analysis.leaflet.LeafletFinder(u22, 'name P*')

bilayer22 = []

for ts in u22.trajectory :
	bilayer22.append(bilayer_thickness(u22, L22))

mean22 = np.mean(bilayer22[25:75])

u23 = MDAnalysis.Universe('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L23 = MDAnalysis.analysis.leaflet.LeafletFinder(u23, 'name P*')

bilayer23 = []

for ts in u23.trajectory :
	bilayer23.append(bilayer_thickness(u23, L23))

mean23 = np.mean(bilayer23[25:75])

av10 = (mean22+mean23)/2

### -60b

u24 = MDAnalysis.Universe('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc')
L24 = MDAnalysis.analysis.leaflet.LeafletFinder(u24, 'name P*')

bilayer24 = []

for ts in u24.trajectory :
	bilayer24.append(bilayer_thickness(u24, L24))

mean24 = np.mean(bilayer24[25:75])

u25 = MDAnalysis.Universe('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc')
L25 = MDAnalysis.analysis.leaflet.LeafletFinder(u25, 'name P*')

bilayer25 = []

for ts in u25.trajectory :
	bilayer25.append(bilayer_thickness(u25, L25))

mean25 = np.mean(bilayer25[25:75])

av11 = (mean24+mean25)/2

##### POPC
### 1b

u26 = MDAnalysis.Universe('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc')
L26 = MDAnalysis.analysis.leaflet.LeafletFinder(u26, 'name P*')

bilayer26 = []

for ts in u26.trajectory :
	bilayer26.append(bilayer_thickness(u26, L26))

mean26 = np.mean(bilayer26[25:75])

u27 = MDAnalysis.Universe('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc')
L27 = MDAnalysis.analysis.leaflet.LeafletFinder(u27, 'name P*')

bilayer27 = []

for ts in u27.trajectory :
	bilayer27.append(bilayer_thickness(u27, L27))

mean27 = np.mean(bilayer27[25:75])

av12 = (mean26+mean27)/2

### -50b

u28 = MDAnalysis.Universe('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L28 = MDAnalysis.analysis.leaflet.LeafletFinder(u28, 'name P*')

bilayer28 = []

for ts in u28.trajectory :
	bilayer28.append(bilayer_thickness(u28, L28))

mean28 = np.mean(bilayer28[25:75])

u29 = MDAnalysis.Universe('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L29 = MDAnalysis.analysis.leaflet.LeafletFinder(u29, 'name P*')

bilayer29 = []

for ts in u29.trajectory :
	bilayer29.append(bilayer_thickness(u29, L29))

mean29 = np.mean(bilayer29[25:75])

av13 = (mean28+mean29)/2

# BILAYER POSITION

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
protCOM = prot.center_of_mass()
calphas = prot.select_atoms('name CA').positions

##########
# PLOTTING
##########

plt.figure(figsize=(10,10))
ax = plt.subplot(111)

plt.plot(calphas[:,0], calphas[:,2], color = 'purple', linewidth = 2)

plt.axhline((Meandif[0] + protCOM[2]) + mean1/2, label = 'F1b1 Average', color = 'blue', linewidth = 1)
plt.axhline((Meandif[0] + protCOM[2]) - mean1/2, color = 'blue', linewidth = 1)
plt.axhline((Meandif[1] + protCOM[2]) + av1/2, label = 'F1b2 Average', color = 'blue', linewidth = 1, alpha = 0.66)
plt.axhline((Meandif[1] + protCOM[2]) - av1/2, color = 'blue', linewidth = 1, alpha = 0.66)
plt.axhline((Meandif[2] + protCOM[2]) + av2/2, label = 'F30b1 Average', color = 'blue', linewidth = 1, linestyle = '--')
plt.axhline((Meandif[2] + protCOM[2]) - av2/2, color = 'blue', linewidth = 1, linestyle = '--')
plt.axhline((Meandif[3] + protCOM[2]) + av2/2, label = 'F30b2 Average', color = 'blue', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline((Meandif[3] + protCOM[2]) + av2/2, color = 'blue', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline((Meandif[4] + protCOM[2]) + av3/2, label = 'F40b1 Average', color = 'blue', linewidth = 1, linestyle = '-.')
plt.axhline((Meandif[4] + protCOM[2]) - av3/2, color = 'blue', linewidth = 1, linestyle = '-.')
plt.axhline((Meandif[5] + protCOM[2]) + av3/2, label = 'F40b2 Average', color = 'blue', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline((Meandif[5] + protCOM[2]) - av3/2, color = 'blue', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline((Meandif[6] + protCOM[2]) + av3/2, label = 'F40b3 Average', color = 'blue', linewidth = 1, linestyle = '-.', alpha = 0.33)
plt.axhline((Meandif[6] + protCOM[2]) - av3/2, color = 'blue', linewidth = 1, linestyle = '-.', alpha = 0.33)

plt.axhline((Meandif[7] + protCOM[2]) + av4/2, label = 'O1b1 Average', color = 'green', linewidth = 1)
plt.axhline((Meandif[7] + protCOM[2]) - av4/2, color = 'green', linewidth = 1)
plt.axhline((Meandif[8] + protCOM[2]) + av4/2, label = 'O1b2 Average', color = 'green', linewidth = 1, alpha = 0.66)
plt.axhline((Meandif[8] + protCOM[2]) - av4/2, color = 'green', linewidth = 1, alpha = 0.66)
plt.axhline((Meandif[9] + protCOM[2]) + av5/2, label = 'O30b1 Average', color = 'green', linewidth = 1, linestyle = '--')
plt.axhline((Meandif[9] + protCOM[2]) - av5/2, color = 'green', linewidth = 1, linestyle = '--')
plt.axhline((Meandif[10] + protCOM[2]) + av5/2, label = 'O30b2 Average', color = 'green', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline((Meandif[10] + protCOM[2]) - av5/2, color = 'green', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline((Meandif[11] + protCOM[2]) + av6/2, label = 'O40b1 Average', color = 'green', linewidth = 1, linestyle = '-.')
plt.axhline((Meandif[11] + protCOM[2]) - av6/2, color = 'green', linewidth = 1, linestyle = '-.')
plt.axhline((Meandif[12] + protCOM[2]) + av6/2, label = 'O40b2 Average', color = 'green', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline((Meandif[12] + protCOM[2]) - av6/2, color = 'green', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline((Meandif[13] + protCOM[2]) + av6/2, label = 'O40b3 Average', color = 'green', linewidth = 1, linestyle = '-.', alpha = 0.33)
plt.axhline((Meandif[13] + protCOM[2]) - av6/2, color = 'green', linewidth = 1, linestyle = '-.', alpha = 0.33)
plt.axhline((Meandif[14] + protCOM[2]) + mean15/2, label = 'O50b1 Average', color = 'green', linewidth = 1, linestyle = ':')
plt.axhline((Meandif[14] + protCOM[2]) - mean15/2, color = 'green', linewidth = 1, linestyle = ':')
plt.axhline((Meandif[15] + protCOM[2]) + av7/2, label = 'O50b2 Average', color = 'green', linewidth = 1, linestyle = ':', alpha = 0.66)
plt.axhline((Meandif[15] + protCOM[2]) - av7/2, color = 'green', linewidth = 1, linestyle = ':', alpha = 0.66)

plt.axhline((Meandif[16] + protCOM[2]) + av8/2, label = 'V1b1 Average', color = 'red', linewidth = 1)
plt.axhline((Meandif[16] + protCOM[2]) - av8/2, color = 'red', linewidth = 1)
plt.axhline((Meandif[17] + protCOM[2]) + av8/2, label = 'V1b2 Average', color = 'red', linewidth = 1, alpha = 0.66)
plt.axhline((Meandif[17] + protCOM[2]) - av8/2, color = 'red', linewidth = 1, alpha = 0.66)
plt.axhline((Meandif[18] + protCOM[2]) + av9/2, label = 'V40b1 Average', color = 'red', linewidth = 1, linestyle = '--')
plt.axhline((Meandif[18] + protCOM[2]) - av9/2, color = 'red', linewidth = 1, linestyle = '--')
plt.axhline((Meandif[19] + protCOM[2]) + av9/2, label = 'V40b2 Average', color = 'red', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline((Meandif[19] + protCOM[2]) - av9/2, color = 'red', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline((Meandif[20] + protCOM[2]) + av9/2, label = 'V40b3 Average', color = 'red', linewidth = 1, linestyle = '--', alpha = 0.33)
plt.axhline((Meandif[20] + protCOM[2]) - av9/2, color = 'red', linewidth = 1, linestyle = '--', alpha = 0.33)
plt.axhline((Meandif[21] + protCOM[2]) + av10/2, label = 'V50b1 Average', color = 'red', linewidth = 1, linestyle = '-.')
plt.axhline((Meandif[21] + protCOM[2]) - av10/2, color = 'red', linewidth = 1, linestyle = '-.')
plt.axhline((Meandif[22] + protCOM[2]) + av10/2, label = 'V50b2 Average', color = 'red', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline((Meandif[22] + protCOM[2]) - av10/2, color = 'red', linewidth = 1, linestyle = '-.', alpha = 0.66)
plt.axhline((Meandif[23] + protCOM[2]) + av11/2, label = 'V60b1 Average', color = 'red', linewidth = 1, linestyle = ':')
plt.axhline((Meandif[23] + protCOM[2]) - av11/2, color = 'red', linewidth = 1, linestyle = ':')
plt.axhline((Meandif[24] + protCOM[2]) + av11/2, label = 'V60b2 Average', color = 'red', linewidth = 1, linestyle = ':', alpha = 0.66)
plt.axhline((Meandif[24] + protCOM[2]) - av11/2, color = 'red', linewidth = 1, linestyle = ':', alpha = 0.66)

plt.axhline((Meandif[25] + protCOM[2]) + av12/2, label = 'P1b1 Average', color = 'black', linewidth = 1)
plt.axhline((Meandif[25] + protCOM[2]) - av12/2, color = 'red', linewidth = 1)
plt.axhline((Meandif[26] + protCOM[2]) + av12/2, label = 'P1b2 Average', color = 'black', linewidth = 1, alpha = 0.66)
plt.axhline((Meandif[26] + protCOM[2]) - av12/2, color = 'red', linewidth = 1, alpha = 0.66)
plt.axhline((Meandif[27] + protCOM[2]) + av13/2, label = 'P50b1 Average', color = 'black', linewidth = 1, linestyle = '--')
plt.axhline((Meandif[27] + protCOM[2]) - av13/2, color = 'red', linewidth = 1, linestyle = '--')
plt.axhline((Meandif[28] + protCOM[2]) + av13/2, label = 'P50b2 Average', color = 'black', linewidth = 1, linestyle = '--', alpha = 0.66)
plt.axhline((Meandif[28] + protCOM[2]) - av13/2, color = 'red', linewidth = 1, linestyle = '--', alpha = 0.66)

plt.ylabel("z coordinate ($\AA$)", fontsize = '18')
plt.xlabel("x coordinate ($\AA$)", fontsize = '18')
w = [20,90,20,120]
plt.axis(w)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(fontsize = 'large', loc = 'center left', bbox_to_anchor=(1,0.5))
plt.title('Bilayer Thickness', fontsize = 30)
plt.savefig("bilayeroverall-prot.png", format='png', dpi=300)
plt.savefig("bilayeroverall-prot.svg", format='svg', dpi=300) 

plt.clf()
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.axhline(av1, label = 'DFPC 1b', color = 'blue', linewidth = 2)
plt.axhline(av2, label = 'DFPC -30b', color = 'blue', linewidth = 2, linestyle = '--')
plt.axhline(av3, label = 'DFPC -40b', color = 'blue', linewidth = 2, linestyle = '-.')
plt.axhline(av4, label = 'DOPC 1b', color = 'green', linewidth = 2)
plt.axhline(av5, label = 'DOPC -30b', color = 'green', linewidth = 2, linestyle = '--')
plt.axhline(av6, label = 'DOPC -40b', color = 'green', linewidth = 2, linestyle = '-.')
plt.axhline(av7, label = 'DOPC -50b', color = 'green', linewidth = 2, linestyle = ':')
plt.axhline(av8, label = 'DVPC 1b', color = 'red', linewidth = 2)
plt.axhline(av9, label = 'DVPC -40b', color = 'red', linewidth = 2, linestyle = '--')
plt.axhline(av10, label = 'DVPC -50b', color = 'red', linewidth = 2, linestyle = '-.')
plt.axhline(av11, label = 'DVPC -60b', color = 'red', linewidth = 2, linestyle = ':')
plt.axhline(av12, label = 'POPC 1b', color = 'black', linewidth = 2)
plt.axhline(av13, label = 'POPC -50b', color = 'black', linewidth = 2, linestyle = '--')

v = [0,1,29,43]
plt.axis(v)
plt.ylabel("Bilayer Thickness ($\AA$)", fontsize = '24')
plt.setp(ax.get_xticklabels(), visible = False)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(fontsize = 'x-large', loc = 'center left', bbox_to_anchor=(1,0.5))
plt.savefig("bilayeroverall.png", format='png', dpi=300)
plt.savefig("bilayeroverall.svg", format='svg', dpi=300)

plt.clf()

dfpc = []
dfpc.append([mean1, mean2, mean3, mean4, mean5, mean6, mean7])

dopc = []
dopc.append([mean8, mean9, mean10, mean11, mean12, mean13, mean14, mean15, mean16])

dvpc = []
dvpc.append([mean17, mean18, mean19, mean20, mean21, mean22, mean23, mean24, mean25])

popc = []
popc.append([mean26, mean27, mean28, mean29])

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
#plt.ylabel('Bilayer Thickness ($\AA$)')
#plt.xlabel('Lipid Species')
#plt.savefig("bilayer-whiskers.png", format='png', dpi=300)
#plt.savefig("bilayer-whiskers.svg", format='svg', dpi=300) 

#plt.clf()
#plt.figure(figsize=(8,8))
#sns.violinplot(data=df, cut = 0)
#plt.ylabel('Bilayer Thickness ($\AA$)')
#plt.xlabel('Lipid Species')
#plt.savefig("bilayer-violin.png", format='png', dpi=300)
#plt.savefig("bilayer-violin.svg", format='svg', dpi=300) 

#plt.clf()
#plt.figure(figsize=(8,8))
#sns.boxplot(data=DF, whis = 100)
#plt.ylabel('Bilayer Thickness ($\AA$)')
#plt.xlabel('Lipid Species')
#plt.savefig("bilayer-whiskers-nopop.png", format='png', dpi=300)
#plt.savefig("bilayer-whiskers-nopop.svg", format='svg', dpi=300) 

#plt.clf()
#plt.figure(figsize=(8,8))
#sns.violinplot(data=DF, cut = 0)
#plt.ylabel('Bilayer Thickness ($\AA$)')
#plt.xlabel('Lipid Species')
#plt.savefig("bilayer-violin-nopop.png", format='png', dpi=300)
#plt.savefig("bilayer-violin-nopop.svg", format='svg', dpi=300) 

plt.figure(figsize=(6,8))
sns.stripplot(data=DF, jitter = 0, size = 10, linewidth = 1)
plt.ylabel('Bilayer Thickness ($\AA$)')
plt.xlabel('Lipid Species')
plt.savefig("bilayer-strip-nopop.png", format='png', dpi=300)
plt.savefig("bilayer-strip-nopop.svg", format='svg', dpi=300)
