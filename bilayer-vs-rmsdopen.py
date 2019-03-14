import numpy as np
import MDAnalysis
import MDAnalysis.analysis.leaflet
import MDAnalysis.analysis.rms
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

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

def bilayer_thickness (A, B) : # finds the bilayer thickness from a file or universe input
	downP = B.groups(0)
	upP = B.groups(1)
	upPxyz = upP.coordinates()
	upPz = upPxyz[:,2] #selects z coordinates
	upPz_av = np.mean(upPz)
	downPxyz = downP.coordinates()
	downPz = downPxyz[:,2] 
	downPz_av = np.mean(downPz)
	bilayer_thickness = downPz_av - upPz_av #difference in average upper and lower z coordinates
	return (bilayer_thickness)

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

########## BILAYER THICKNESS

thickness = []

##### DFPC
### 1b

u1 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc')
L1 = MDAnalysis.analysis.leaflet.LeafletFinder(u1, 'name P*')

bilayer1 = []

for ts in u1.trajectory :
	bilayer1.append(bilayer_thickness(u1, L1))

thickness.append(np.mean(bilayer1[25:75]))

#u2 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
#L2 = MDAnalysis.analysis.leaflet.LeafletFinder(u2, 'name P*')

#bilayer2 = []

#for ts in u2.trajectory :
#	bilayer2.append(bilayer_thickness(u2, L2))

#thickness.append(np.mean(bilayer2[25:75]))

### -30b

u3 = MDAnalysis.Universe('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc')
L3 = MDAnalysis.analysis.leaflet.LeafletFinder(u3, 'name P*')

bilayer3 = []

for ts in u3.trajectory :
	bilayer3.append(bilayer_thickness(u3, L3))

thickness.append(np.mean(bilayer3[25:75]))

u4 = MDAnalysis.Universe('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc')
L4 = MDAnalysis.analysis.leaflet.LeafletFinder(u4, 'name P*')

bilayer4 = []

for ts in u4.trajectory :
	bilayer4.append(bilayer_thickness(u4, L4))

thickness.append(np.mean(bilayer4[25:75]))

### -40b

bilayer5 = []

u5 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L5 = MDAnalysis.analysis.leaflet.LeafletFinder(u5, 'name P*') #selects P atoms and puts in two leaflet

for ts in u5.trajectory :
	bilayer5.append(bilayer_thickness(u5, L5))

thickness.append(np.mean(bilayer5[25:75]))

bilayer6 = []

u6 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L6 = MDAnalysis.analysis.leaflet.LeafletFinder(u6, 'name P*') #selects P atoms and puts in two leaflet

for ts in u6.trajectory :
	bilayer6.append(bilayer_thickness(u6, L6))

thickness.append(np.mean(bilayer6[25:75]))

bilayer7 = []

u7 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L7 = MDAnalysis.analysis.leaflet.LeafletFinder(u7, 'name P*') #selects P atoms and puts in two leaflet

for ts in u7.trajectory :
	bilayer7.append(bilayer_thickness(u7, L7))

thickness.append(np.mean(bilayer7[25:75]))

##### DOPC
### 1b

u8 = MDAnalysis.Universe('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc')
L8 = MDAnalysis.analysis.leaflet.LeafletFinder(u8, 'name P*')

bilayer8 = []

for ts in u8.trajectory :
	bilayer8.append(bilayer_thickness(u8, L8))

thickness.append(np.mean(bilayer8[25:75]))

u9 = MDAnalysis.Universe('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L9 = MDAnalysis.analysis.leaflet.LeafletFinder(u9, 'name P*')

bilayer9 = []

for ts in u9.trajectory :
	bilayer9.append(bilayer_thickness(u9, L9))

thickness.append(np.mean(bilayer9[25:75]))

### -30b

u10 = MDAnalysis.Universe('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc')
L10 = MDAnalysis.analysis.leaflet.LeafletFinder(u10, 'name P*')

bilayer10 = []

for ts in u10.trajectory :
	bilayer10.append(bilayer_thickness(u10, L10))

thickness.append(np.mean(bilayer10[25:75]))

u11 = MDAnalysis.Universe('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc')
L11 = MDAnalysis.analysis.leaflet.LeafletFinder(u11, 'name P*')

bilayer11 = []

for ts in u11.trajectory :
	bilayer11.append(bilayer_thickness(u11, L11))

thickness.append(np.mean(bilayer11[25:75]))

### -40b

bilayer12 = []

u12 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L12 = MDAnalysis.analysis.leaflet.LeafletFinder(u12, 'name P*') #selects P atoms and puts in two leaflet

for ts in u12.trajectory :
	bilayer12.append(bilayer_thickness(u12, L12))

thickness.append(np.mean(bilayer12[25:75]))

bilayer13 = []

u13 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L13 = MDAnalysis.analysis.leaflet.LeafletFinder(u13, 'name P*') #selects P atoms and puts in two leaflet

for ts in u13.trajectory :
	bilayer13.append(bilayer_thickness(u13, L13))

thickness.append(np.mean(bilayer13[25:75]))

bilayer14 = []

u14 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L14 = MDAnalysis.analysis.leaflet.LeafletFinder(u14, 'name P*') #selects P atoms and puts in two leaflet

for ts in u14.trajectory :
	bilayer14.append(bilayer_thickness(u14, L14))

thickness.append(np.mean(bilayer14[25:75]))

### -50b

u15 = MDAnalysis.Universe('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L15 = MDAnalysis.analysis.leaflet.LeafletFinder(u15, 'name P*')

bilayer15 = []

for ts in u15.trajectory :
	bilayer15.append(bilayer_thickness(u15, L15))

thickness.append(np.mean(bilayer15[25:75]))

u16 = MDAnalysis.Universe('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L16 = MDAnalysis.analysis.leaflet.LeafletFinder(u16, 'name P*')

bilayer16 = []

for ts in u16.trajectory :
	bilayer16.append(bilayer_thickness(u16, L16))

thickness.append(np.mean(bilayer16[25:75]))

##### DVPC
### 1b

u17 = MDAnalysis.Universe('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc')
L17 = MDAnalysis.analysis.leaflet.LeafletFinder(u17, 'name P*')

bilayer17 = []

for ts in u17.trajectory :
	bilayer17.append(bilayer_thickness(u17, L17))

thickness.append(np.mean(bilayer17[25:75]))

u18 = MDAnalysis.Universe('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L18 = MDAnalysis.analysis.leaflet.LeafletFinder(u18, 'name P*')

bilayer18 = []

for ts in u18.trajectory :
	bilayer18.append(bilayer_thickness(u18, L18))

thickness.append(np.mean(bilayer18[25:75]))

### -40b

bilayer19 = []

u19 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L19 = MDAnalysis.analysis.leaflet.LeafletFinder(u19, 'name P*') #selects P atoms and puts in two leaflet

for ts in u19.trajectory :
	bilayer19.append(bilayer_thickness(u19, L19))

thickness.append(np.mean(bilayer19[25:75]))

bilayer20 = []

u20 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L20 = MDAnalysis.analysis.leaflet.LeafletFinder(u20, 'name P*') #selects P atoms and puts in two leaflet

for ts in u20.trajectory :
	bilayer20.append(bilayer_thickness(u20, L20))

thickness.append(np.mean(bilayer20[25:75]))

bilayer21 = []

u21 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L21 = MDAnalysis.analysis.leaflet.LeafletFinder(u21, 'name P*') #selects P atoms and puts in two leaflet

for ts in u21.trajectory :
	bilayer21.append(bilayer_thickness(u21, L21))

thickness.append(np.mean(bilayer21[25:75]))

### -50b

u22 = MDAnalysis.Universe('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L22 = MDAnalysis.analysis.leaflet.LeafletFinder(u22, 'name P*')

bilayer22 = []

for ts in u22.trajectory :
	bilayer22.append(bilayer_thickness(u22, L22))

thickness.append(np.mean(bilayer22[25:75]))

u23 = MDAnalysis.Universe('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L23 = MDAnalysis.analysis.leaflet.LeafletFinder(u23, 'name P*')

bilayer23 = []

for ts in u23.trajectory :
	bilayer23.append(bilayer_thickness(u23, L23))

thickness.append(np.mean(bilayer23[25:75]))

### -60b

u24 = MDAnalysis.Universe('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc')
L24 = MDAnalysis.analysis.leaflet.LeafletFinder(u24, 'name P*')

bilayer24 = []

for ts in u24.trajectory :
	bilayer24.append(bilayer_thickness(u24, L24))

thickness.append(np.mean(bilayer24[25:75]))

u25 = MDAnalysis.Universe('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc')
L25 = MDAnalysis.analysis.leaflet.LeafletFinder(u25, 'name P*')

bilayer25 = []

for ts in u25.trajectory :
	bilayer25.append(bilayer_thickness(u25, L25))

thickness.append(np.mean(bilayer25[25:75]))

########## RMSD

rmsd = []

##### DFPC
### 1b

F1b1RMS = []
F1b1TS = []
RMSD_traj('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc', F1b1RMS, F1b1TS)
rmsd.append(F1b1TS[0])

F1b2RMS = []
F1b2TS = []
RMSD_traj('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', F1b2RMS, F1b2TS)
rmsd.append(F1b2TS[0])

### -30b

F30b1RMS = []
F30b1TS = []
RMSD_traj('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', F30b1RMS, F30b1TS)
rmsd.append(F30b1TS[0])

F30b2RMS = []
F30b2TS = []
RMSD_traj('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', F30b2RMS, F30b2TS)
rmsd.append(F30b2TS[0])

### -40b

F40b1RMS = []
F40b1TS = []
RMSD_traj('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', F40b1RMS, F40b1TS)
rmsd.append(F40b1TS[0])

F40b2RMS = []
F40b2TS = []
RMSD_traj('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', F40b2RMS, F40b2TS)
rmsd.append(F40b2TS[0])

F40b3RMS = []
F40b3TS = []
RMSD_traj('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', F40b3RMS, F40b3TS)
rmsd.append(F40b3TS[0])

##### DOPC
### 1b

O1b1RMS = []
O1b1TS = []
RMSD_traj('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc', O1b1RMS, O1b1TS)
rmsd.append(O1b1TS[0])

O1b2RMS = []
O1b2TS = []
RMSD_traj('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', O1b1RMS, O1b2TS)
rmsd.append(O1b2TS[0])

### -30b

O30b1RMS = []
O30b1TS = []
RMSD_traj('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', O30b1RMS, O30b1TS)
rmsd.append(O30b1TS[0])

O30b2RMS = []
O30b2TS = []
RMSD_traj('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', O30b2RMS, O30b2TS)
rmsd.append(O30b2TS[0])

### -40b

O40b1RMS = []
O40b1TS = []
RMSD_traj('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', O40b1RMS, O40b1TS)
rmsd.append(O40b1TS[0])

O40b2RMS = []
O40b2TS = []
RMSD_traj('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', O40b2RMS, O40b2TS)
rmsd.append(O40b2TS[0])

O40b3RMS = []
O40b3TS = []
RMSD_traj('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', O40b3RMS, O40b3TS)
rmsd.append(O40b3TS[0])

### -50b

O50b1RMS = []
O50b1TS = []
RMSD_traj('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', O50b1RMS, O50b1TS)
rmsd.append(O50b1TS[0])

O50b2RMS = []
O50b2TS = []
RMSD_traj('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', O50b2RMS, O50b2TS)
rmsd.append(O50b2TS[0])

##### DVPC
### 1b

V1b1RMS = []
V1b1TS = []
RMSD_traj('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc', V1b1RMS, V1b1TS)
rmsd.append(V1b1TS[0])

V1b2RMS = []
V1b2TS = []
RMSD_traj('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', V1b2RMS, V1b2TS)
rmsd.append(V1b2TS[0])

### -40b

V40b1RMS = []
V40b1TS = []
RMSD_traj('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', V40b1RMS, V40b1TS)
rmsd.append(V40b1TS[0])

V40b2RMS = []
V40b2TS = []
RMSD_traj('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', V40b2RMS, V40b2TS)
rmsd.append(V40b2TS[0])

V40b3RMS = []
V40b3TS = []
RMSD_traj('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', V40b3RMS, V40b3TS)
rmsd.append(V40b3TS[0])

### -50b

V50b1RMS = []
V50b1TS = []
RMSD_traj('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', V50b1RMS, V50b1TS)
rmsd.append(V50b1TS[0])

V50b2RMS = []
V50b2TS = []
RMSD_traj('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', V50b2RMS, V50b2TS)
rmsd.append(V50b2TS[0])

### -60b

V60b1RMS = []
V60b1TS = []
RMSD_traj('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', V60b1RMS, V60b1TS)
rmsd.append(V60b1TS[0])

V60b2RMS = []
V60b2TS = []
RMSD_traj('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', V60b2RMS, V60b2TS)
rmsd.append(V60b2TS[0]) 

##########
# PLOTTING
##########

plt.figure(figsize = (12,8))
plt.scatter(rmsd, thickness)
v = [0,100,29,36]
plt.axis(v)
plt.ylabel("Bilayer Thickness ($\AA$)", fontsize = '24')
plt.xlabel("Opening Time (ns)", fontsize = '24')
plt.savefig("bilayervsrmsd.png", format='png', dpi=300)
plt.savefig("bilayervsrmsd.svg", format='svg', dpi=300)                               
