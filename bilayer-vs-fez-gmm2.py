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
# DEFINITIONS
##########

def dis (sel1, sel2, x) :
	res_i = x.select_atoms(sel1 + ' and name CA').positions
	res_j = x.select_atoms(sel2 + ' and name CA').positions
	dis = np.linalg.norm(res_i - res_j)
	return dis

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
		if i < 6.1 :
			Fen_A_2.append(Fen_A.index(i))
	for i in Fen_B :
		if i < 6.1 :
			Fen_B_2.append(Fen_B.index(i))
	for i in Zip_A :
		if i > 11.2 :
			Zip_A_2.append(Zip_A.index(i))
	for i in Zip_B :
		if i > 11.2 :
			Zip_B_2.append(Zip_B.index(i))
	for i in Exp_A :
		if i > 11.2 :
			Exp_A_2.append(Exp_A.index(i))
	for i in Exp_B :
		if i > 11.2 :
			Exp_B_2.append(Exp_B.index(i))
	for i in Fen_A_2 :
		if i in Zip_A_2 and Exp_A_2 :
			C.append(i)
	for i in Fen_B_2 :
		if i in Zip_B_2 and Exp_B_2 :
			D.append(i)
	if C == [] :
		C.append(float(101))
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

u2 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L2 = MDAnalysis.analysis.leaflet.LeafletFinder(u2, 'name P*')

bilayer2 = []

for ts in u2.trajectory :
	bilayer2.append(bilayer_thickness(u2, L2))

thickness.append(np.mean(bilayer2[25:75]))

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

##### POPC
### 1b

u26 = MDAnalysis.Universe('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc')
L26 = MDAnalysis.analysis.leaflet.LeafletFinder(u26, 'name P*')

bilayer26 = []

for ts in u26.trajectory :
	bilayer26.append(bilayer_thickness(u26, L26))

thickness.append(np.mean(bilayer26[25:75]))

u27 = MDAnalysis.Universe('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc')
L27 = MDAnalysis.analysis.leaflet.LeafletFinder(u27, 'name P*')

bilayer27 = []

for ts in u27.trajectory :
	bilayer27.append(bilayer_thickness(u27, L27))

thickness.append(np.mean(bilayer27[25:75]))

### -50b

u28 = MDAnalysis.Universe('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L28 = MDAnalysis.analysis.leaflet.LeafletFinder(u28, 'name P*')

bilayer28 = []

for ts in u28.trajectory :
	bilayer28.append(bilayer_thickness(u28, L28))

thickness.append(np.mean(bilayer28[25:75]))

u29 = MDAnalysis.Universe('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L29 = MDAnalysis.analysis.leaflet.LeafletFinder(u29, 'name P*')

bilayer29 = []

for ts in u29.trajectory :
	bilayer29.append(bilayer_thickness(u29, L29))

thickness.append(np.mean(bilayer29[25:75]))

########## FEZ DISTANCES

fezA = []
fezB = []

##### DFPC
### 1b

F1b1chA = []
F1b1chB = []
FEZ_traj('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc', F1b1chA, F1b1chB)
fezA.append(F1b1chA[0])
fezB.append(F1b1chB[0])

F1b2chA = []
F1b2chB = []
FEZ_traj('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', F1b2chA, F1b2chB)
fezA.append(F1b2chA[0])
fezB.append(F1b2chB[0])

### -30b

F30b1chA = []
F30b1chB = []
FEZ_traj('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', F30b1chA, F30b1chB)
fezA.append(F30b1chA[0])
fezB.append(F30b1chB[0])

F30b2chA = []
F30b2chB = []
FEZ_traj('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', F30b2chA, F30b2chB)
fezA.append(F30b2chA[0])
fezB.append(F30b2chB[0])

### -40b

F40b1chA = []
F40b1chB = []
FEZ_traj('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', F40b1chA, F40b1chB)
fezA.append(F40b1chA[0])
fezB.append(F40b1chB[0])

F40b2chA = []
F40b2chB = []
FEZ_traj('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', F40b2chA, F40b2chB)
fezA.append(F40b2chA[0])
fezB.append(F40b2chB[0])

F40b3chA = []
F40b3chB = []
FEZ_traj('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', F40b3chA, F40b3chB)
fezA.append(F40b3chA[0])
fezB.append(F40b3chB[0])

##### DOPC
### 1b

O1b1chA = []
O1b1chB = []
FEZ_traj('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc', O1b1chA, O1b1chB)
fezA.append(O1b1chA[0])
fezB.append(O1b1chB[0])

O1b2chA = []
O1b2chB = []
FEZ_traj('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', O1b2chA, O1b2chB)
fezA.append(O1b2chA[0])
fezB.append(O1b2chB[0])

### -30b

O30b1chA = []
O30b1chB = []
FEZ_traj('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', O30b1chA, O30b1chB)
fezA.append(O30b1chA[0])
fezB.append(O30b1chB[0])

O30b2chA = []
O30b2chB = []
FEZ_traj('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', O30b2chA, O30b2chB)
fezA.append(O30b2chA[0])
fezB.append(O30b2chB[0])

### -40b

O40b1chA = []
O40b1chB = []
FEZ_traj('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', O40b1chA, O40b1chB)
fezA.append(O40b1chA[0])
fezB.append(O40b1chB[0])

O40b2chA = []
O40b2chB = []
FEZ_traj('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', O40b2chA, O40b2chB)
fezA.append(O40b2chA[0])
fezB.append(O40b2chB[0])

O40b3chA = []
O40b3chB = []
FEZ_traj('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', O40b3chA, O40b3chB)
fezA.append(O40b3chA[0])
fezB.append(O40b3chB[0])

### -50b

O50b1chA = []
O50b1chB = []
FEZ_traj('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', O50b1chA, O50b1chB)
fezA.append(O50b1chA[0])
fezB.append(O50b1chB[0])

O50b2chA = []
O50b2chB = []
FEZ_traj('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', O50b2chA, O50b2chB)
fezA.append(O50b2chA[0])
fezB.append(O50b2chB[0])

##### DVPC
### 1b

V1b1chA = []
V1b1chB = []
FEZ_traj('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc', V1b1chA, V1b1chB)
fezA.append(V1b1chA[0])
fezB.append(V1b1chB[0])

V1b2chA = []
V1b2chB = []
FEZ_traj('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', V1b2chA, V1b2chB)
fezA.append(V1b2chA[0])
fezB.append(V1b2chB[0])

### -40b

V40b1chA = []
V40b1chB = []
FEZ_traj('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', V40b1chA, V40b1chB)
fezA.append(V40b1chA[0])
fezB.append(V40b1chB[0])

V40b2chA = []
V40b2chB = []
FEZ_traj('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', V40b2chA, V40b2chB)
fezA.append(V40b2chA[0])
fezB.append(V40b2chB[0])

V40b3chA = []
V40b3chB = []
FEZ_traj('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', V40b3chA, V40b3chB)
fezA.append(V40b3chA[0])
fezB.append(V40b3chB[0])

### -50b

V50b1chA = []
V50b1chB = []
FEZ_traj('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', V50b1chA, V50b1chB)
fezA.append(V50b1chA[0])
fezB.append(V50b1chB[0])

V50b2chA = []
V50b2chB = []
FEZ_traj('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', V50b2chA, V50b2chB)
fezA.append(V50b2chA[0])
fezB.append(V50b2chB[0])

### -60b

V60b1chA = []
V60b1chB = []
FEZ_traj('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', V60b1chA, V60b1chB)
fezA.append(V60b1chA[0])
fezB.append(V60b1chB[0])

V60b2chA = []
V60b2chB = []
FEZ_traj('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', V60b2chA, V60b2chB)
fezA.append(V60b2chA[0])
fezB.append(V60b2chB[0])

##### POPC
### 1b

P1b1chA = []
P1b1chB = []
FEZ_traj('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc', P1b1chA, P1b1chB)
fezA.append(P1b1chA[0])
fezB.append(P1b1chB[0])

P1b2chA = []
P1b2chB = []
FEZ_traj('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc', P1b2chA, P1b2chB)
fezA.append(P1b2chA[0])
fezB.append(P1b2chB[0])

### -50b

P50b1chA = []
P50b1chB = []
FEZ_traj('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', P50b1chA, P50b1chB)
fezA.append(P50b1chA[0])
fezB.append(P50b1chB[0])

P50b2chA = []
P50b2chB = []
FEZ_traj('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', P50b2chA, P50b2chB)
fezA.append(P50b2chA[0])
fezB.append(P50b2chB[0])

##########
# PLOTTING
##########

fig = plt.figure(figsize=(16,8))
ax = plt.subplot(111)
plt.scatter(fezA[0], thickness[0], label = 'DFPC 1b Run 1', color = '#9647ef', marker = '^')
plt.scatter(fezB[0], thickness[0], color = '#9647ef', marker = '^')
plt.scatter(fezA[1], thickness[1], label = 'DFPC 1b Run 2', color = '#9647ef', marker = '^', alpha = 0.66)
plt.scatter(fezB[1], thickness[1], color = '#9647ef', marker = '^', alpha = 0.66)
plt.scatter(fezA[2], thickness[2], label = 'DFPC -30b Run 1', color = '#9647ef', marker = '<')
plt.scatter(fezB[2], thickness[2], color = '#9647ef', marker = '<')
plt.scatter(fezA[3], thickness[3], label = 'DFPC -30b Run 2', color = '#9647ef', marker = '<', alpha = 0.66)
plt.scatter(fezB[3], thickness[3], color = '#9647ef', marker = '<', alpha = 0.66)
plt.scatter(fezA[4], thickness[4], label = 'DFPC -40b Run 1', color = '#9647ef', marker = '>')
plt.scatter(fezB[4], thickness[4], color = '#9647ef', marker = '>')
plt.scatter(fezA[5], thickness[5], label = 'DFPC -40b Run 2', color = '#9647ef', marker = '>', alpha = 0.66)
plt.scatter(fezB[5], thickness[5], color = '#9647ef', marker = '>', alpha = 0.66)
plt.scatter(fezA[6], thickness[6], label = 'DFPC -40b Run 3', color = '#9647ef', marker = '>', alpha = 0.33)
plt.scatter(fezB[6], thickness[6], color = '#9647ef', marker = '>', alpha = 0.33)
plt.scatter(fezA[7], thickness[7], label = 'DOPC 1b Run 1', color = '#10c8aa', marker = '^')
plt.scatter(fezB[7], thickness[7], color = '#10c8aa', marker = '^')
plt.scatter(fezA[8], thickness[8], label = 'DOPC 1b Run 2', color = '#10c8aa', marker = '^', alpha = 0.66)
plt.scatter(fezB[8], thickness[8], color = '#10c8aa', marker = '^', alpha = 0.66)
plt.scatter(fezA[9], thickness[9], label = 'DOPC -30b Run 1', color = '#10c8aa', marker = '<')
plt.scatter(fezB[9], thickness[9], color = '#10c8aa', marker = '<')
plt.scatter(fezA[10], thickness[10], label = 'DOPC -30b Run 2', color = '#10c8aa', marker = '<', alpha = 0.66)
plt.scatter(fezB[10], thickness[10], color = '#10c8aa', marker = '<', alpha = 0.66)
plt.scatter(fezA[11], thickness[11], label = 'DOPC -40b Run 1', color = '#10c8aa', marker = '>')
plt.scatter(fezB[11], thickness[11], color = '#10c8aa', marker = '>')
plt.scatter(fezA[12], thickness[12], label = 'DOPC -40b Run 2', color = '#10c8aa', marker = '>', alpha = 0.66)
plt.scatter(fezB[12], thickness[12], color = '#10c8aa', marker = '>', alpha = 0.66)
plt.scatter(fezA[13], thickness[13], label = 'DOPC -40b Run 3', color = '#10c8aa', marker = '>', alpha = 0.33)
plt.scatter(fezB[13], thickness[13], color = '#10c8aa', marker = '>', alpha = 0.33)
plt.scatter(fezA[14], thickness[14], label = 'DOPC -50b Run 1', color = '#10c8aa', marker = 'v')
plt.scatter(fezB[14], thickness[14], color = '#10c8aa', marker = 'v')
plt.scatter(fezA[15], thickness[15], label = 'DOPC -50b Run 2', color = '#10c8aa', marker = 'v', alpha = 0.66)
plt.scatter(fezB[15], thickness[15], color = '#10c8aa', marker = 'v', alpha = 0.66)
plt.scatter(fezA[16], thickness[16], label = 'DVPC 1b Run 1', color = '#f9c116', marker = '^')
plt.scatter(fezB[16], thickness[16], color = '#f9c116', marker = '^')
plt.scatter(fezA[17], thickness[17], label = 'DVPC 1b Run 2', color = '#f9c116', marker = '^', alpha = 0.66)
plt.scatter(fezB[17], thickness[17], color = '#f9c116', marker = '^', alpha = 0.66)
plt.scatter(fezA[18], thickness[18], label = 'DVPC -40b Run 1', color = '#f9c116', marker = '<')
plt.scatter(fezB[18], thickness[18], color = '#f9c116', marker = '<')
plt.scatter(fezA[19], thickness[19], label = 'DVPC -40b Run 2', color = '#f9c116', marker = '<', alpha = 0.66)
plt.scatter(fezB[19], thickness[19], color = '#f9c116', marker = '<', alpha = 0.66)
plt.scatter(fezA[20], thickness[20], label = 'DVPC -40b Run 3', color = '#f9c116', marker = '<', alpha = 0.33)
plt.scatter(fezB[20], thickness[20], color = '#f9c116', marker = '<', alpha = 0.33)
plt.scatter(fezA[21], thickness[21], label = 'DVPC -50b Run 1', color = '#f9c116', marker = '>')
plt.scatter(fezB[21], thickness[21], color = '#f9c116', marker = '>')
plt.scatter(fezA[22], thickness[22], label = 'DVPC -50b Run 2', color = '#f9c116', marker = '>', alpha = 0.66)
plt.scatter(fezB[22], thickness[22], color = '#f9c116', marker = '>', alpha = 0.66)
plt.scatter(fezA[23], thickness[23], label = 'DVPC -60b Run 1', color = '#f9c116', marker = 'v')
plt.scatter(fezB[23], thickness[23], color = '#f9c116', marker = 'v')
plt.scatter(fezA[24], thickness[24], label = 'DVPC -60b Run 2', color = '#f9c116', marker = 'v', alpha = 0.66)
plt.scatter(fezB[24], thickness[24], color = '#f9c116', marker = 'v', alpha = 0.66)
#plt.scatter(fezA[25], thickness[25], label = 'POPC 1b Run 1', color = 'black', marker = '^')
#plt.scatter(fezB[25], thickness[25], color = 'black', marker = '^')
#plt.scatter(fezA[26], thickness[26], label = 'POPC 1b Run 2', color = 'black', marker = '^', alpha = 0.66)
#plt.scatter(fezB[26], thickness[26], color = 'black', marker = '^', alpha = 0.66)
#plt.scatter(fezA[27], thickness[27], label = 'POPC -50b Run 1', color = 'black', marker = '<')
#plt.scatter(fezB[27], thickness[27], color = 'black', marker = '<')
#plt.scatter(fezA[28], thickness[28], label = 'POPC -50b Run 2', color = 'black', marker = '<', alpha = 0.66)
#plt.scatter(fezB[28], thickness[28], color = 'black', marker = '<', alpha = 0.66)

v = [0,102,29,36]
plt.axis(v)
plt.axvspan(100, 102, facecolor = 'black', alpha = 0.2)
plt.ylabel("Bilayer Thickness ($\AA$)", fontsize = '24')
plt.xlabel("Opening Time (ns)", fontsize = '24')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(fontsize = 'small', loc = 'center left', bbox_to_anchor=(1.05,0.5))
plt.title("Bilayer Thickness vs Time when F < 6.1, E > 11.2, Z > 11.2 ($\AA$)", fontsize = 24)
plt.savefig("bilayervsfez-gmm2.png", format='png', dpi=300)
plt.savefig("bilayervsfez-gmm2.svg", format='svg', dpi=300)                               
