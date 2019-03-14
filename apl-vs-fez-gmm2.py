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

########## APL

apl = []

##### DFPC
### 1b

d1 = np.loadtxt('../../DFPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)
mean1 = np.mean(d1[25:75, 1])
apl.append(mean1)

d2 = np.loadtxt('../../DFPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)
mean2 = np.mean(d2[25:75, 1])
apl.append(mean2)

### -30b

d3 = np.loadtxt('../../DFPC/Backwards/md/30bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean3 = np.mean(d3[25:75, 1])
apl.append(mean3)

d4 = np.loadtxt('../../DFPC/Backwards/md/30bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean4 = np.mean(d4[25:75, 1])
apl.append(mean4)

### -40b

d5 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean5 = np.mean(d5[25:75, 1])
apl.append(mean5)

d6 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean6 = np.mean(d6[25:75, 1])
apl.append(mean6)

d7 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)
mean7 = np.mean(d7[25:75, 1])
apl.append(mean7)

##### DOPC
### 1b

d8 = np.loadtxt('../../DOPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)
mean8 = np.mean(d8[25:75, 1])
apl.append(mean8)

d9 = np.loadtxt('../../DOPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)
mean9 = np.mean(d9[25:75, 1])
apl.append(mean9)

### -30b

d10 = np.loadtxt('../../DOPC/Backwards/md/30bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean10 = np.mean(d10[25:75, 1])
apl.append(mean10)

d11 = np.loadtxt('../../DOPC/Backwards/md/30bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean11 = np.mean(d11[25:75, 1])
apl.append(mean11)

### -40b

d12 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean12 = np.mean(d12[25:75, 1])
apl.append(mean12)

d13 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean13 = np.mean(d13[25:75, 1])
apl.append(mean13)

d14 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)
mean14 = np.mean(d14[25:75, 1])
apl.append(mean14)

### -50b

d15 = np.loadtxt('../../DOPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean15 = np.mean(d15[25:75, 1])
apl.append(mean15)

d16 = np.loadtxt('../../DOPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean16 = np.mean(d16[25:75, 1])
apl.append(mean16)

##### DVPC
### 1b

d17 = np.loadtxt('../../DVPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)
mean17 = np.mean(d17[25:75, 1])
apl.append(mean17)

d18 = np.loadtxt('../../DVPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)
mean18 = np.mean(d18[25:75, 1])
apl.append(mean18)

### -40b

d19 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean19 = np.mean(d19[25:75, 1])
apl.append(mean19)

d20 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean20 = np.mean(d20[25:75, 1])
apl.append(mean20)

d21 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)
mean21 = np.mean(d21[25:75, 1])
apl.append(mean21)

### -50b

d22 = np.loadtxt('../../DVPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean22 = np.mean(d22[25:75, 1])
apl.append(mean22)

d23 = np.loadtxt('../../DVPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean23 = np.mean(d23[25:75, 1])
apl.append(mean23)

### -60b

d24 = np.loadtxt('../../DVPC/Backwards/md/60bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean24 = np.mean(d24[25:75, 1])
apl.append(mean24)

d25 = np.loadtxt('../../DVPC/Backwards/md/60bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean25 = np.mean(d25[25:75, 1])
apl.append(mean25)

##### POPC
### 1b

d26 = np.loadtxt('../../POPC/Backwards/md/1bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean26 = np.mean(d26[25:75, 1])
apl.append(mean26)

d27 = np.loadtxt('../../POPC/Backwards/md/1bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean27 = np.mean(d27[25:75, 1])
apl.append(mean27)

### -50b

d28 = np.loadtxt('../../POPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)
mean28 = np.mean(d28[25:75, 1])
apl.append(mean28)

d29 = np.loadtxt('../../POPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)
mean29 = np.mean(d29[25:75, 1])
apl.append(mean29)

###

apl = apl*100

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
plt.scatter(fezA[0], apl[0]*100, label = 'DFPC 1b Run 1', color = '#9647ef', marker = '^')
plt.scatter(fezB[0], apl[0]*100, color = '#9647ef', marker = '^')
plt.scatter(fezA[1], apl[1]*100, label = 'DFPC 1b Run 2', color = '#9647ef', marker = '^', alpha = 0.66)
plt.scatter(fezB[1], apl[1]*100, color = '#9647ef', marker = '^', alpha = 0.66)
plt.scatter(fezA[2], apl[2]*100, label = 'DFPC -30b Run 1', color = '#9647ef', marker = '<')
plt.scatter(fezB[2], apl[2]*100, color = '#9647ef', marker = '<')
plt.scatter(fezA[3], apl[3]*100, label = 'DFPC -30b Run 2', color = '#9647ef', marker = '<', alpha = 0.66)
plt.scatter(fezB[3], apl[3]*100, color = '#9647ef', marker = '<', alpha = 0.66)
plt.scatter(fezA[4], apl[4]*100, label = 'DFPC -40b Run 1', color = '#9647ef', marker = '>')
plt.scatter(fezB[4], apl[4]*100, color = '#9647ef', marker = '>')
plt.scatter(fezA[5], apl[5]*100, label = 'DFPC -40b Run 2', color = '#9647ef', marker = '>', alpha = 0.66)
plt.scatter(fezB[5], apl[5]*100, color = '#9647ef', marker = '>', alpha = 0.66)
plt.scatter(fezA[6], apl[6]*100, label = 'DFPC -40b Run 3', color = '#9647ef', marker = '>', alpha = 0.33)
plt.scatter(fezB[6], apl[6]*100, color = '#9647ef', marker = '>', alpha = 0.33)
plt.scatter(fezA[7], apl[7]*100, label = 'DOPC 1b Run 1', color = '#10c8aa', marker = '^')
plt.scatter(fezB[7], apl[7]*100, color = '#10c8aa', marker = '^')
plt.scatter(fezA[8], apl[8]*100, label = 'DOPC 1b Run 2', color = '#10c8aa', marker = '^', alpha = 0.66)
plt.scatter(fezB[8], apl[8]*100, color = '#10c8aa', marker = '^', alpha = 0.66)
plt.scatter(fezA[9], apl[9]*100, label = 'DOPC -30b Run 1', color = '#10c8aa', marker = '<')
plt.scatter(fezB[9], apl[9]*100, color = '#10c8aa', marker = '<')
plt.scatter(fezA[10], apl[10]*100, label = 'DOPC -30b Run 2', color = '#10c8aa', marker = '<', alpha = 0.66)
plt.scatter(fezB[10], apl[10]*100, color = '#10c8aa', marker = '<', alpha = 0.66)
plt.scatter(fezA[11], apl[11]*100, label = 'DOPC -40b Run 1', color = '#10c8aa', marker = '>')
plt.scatter(fezB[11], apl[11]*100, color = '#10c8aa', marker = '>')
plt.scatter(fezA[12], apl[12]*100, label = 'DOPC -40b Run 2', color = '#10c8aa', marker = '>', alpha = 0.66)
plt.scatter(fezB[12], apl[12]*100, color = '#10c8aa', marker = '>', alpha = 0.66)
plt.scatter(fezA[13], apl[13]*100, label = 'DOPC -40b Run 3', color = '#10c8aa', marker = '>', alpha = 0.33)
plt.scatter(fezB[13], apl[13]*100, color = '#10c8aa', marker = '>', alpha = 0.33)
plt.scatter(fezA[14], apl[14]*100, label = 'DOPC -50b Run 1', color = '#10c8aa', marker = 'v')
plt.scatter(fezB[14], apl[14]*100, color = '#10c8aa', marker = 'v')
plt.scatter(fezA[15], apl[15]*100, label = 'DOPC -50b Run 2', color = '#10c8aa', marker = 'v', alpha = 0.66)
plt.scatter(fezB[15], apl[15]*100, color = '#10c8aa', marker = 'v', alpha = 0.66)
plt.scatter(fezA[16], apl[16]*100, label = 'DVPC 1b Run 1', color = '#f9c116', marker = '^')
plt.scatter(fezB[16], apl[16]*100, color = '#f9c116', marker = '^')
plt.scatter(fezA[17], apl[17]*100, label = 'DVPC 1b Run 2', color = '#f9c116', marker = '^', alpha = 0.66)
plt.scatter(fezB[17], apl[17]*100, color = '#f9c116', marker = '^', alpha = 0.66)
plt.scatter(fezA[18], apl[18]*100, label = 'DVPC -40b Run 1', color = '#f9c116', marker = '<')
plt.scatter(fezB[18], apl[18]*100, color = '#f9c116', marker = '<')
plt.scatter(fezA[19], apl[19]*100, label = 'DVPC -40b Run 2', color = '#f9c116', marker = '<', alpha = 0.66)
plt.scatter(fezB[19], apl[19]*100, color = '#f9c116', marker = '<', alpha = 0.66)
plt.scatter(fezA[20], apl[20]*100, label = 'DVPC -40b Run 3', color = '#f9c116', marker = '<', alpha = 0.33)
plt.scatter(fezB[20], apl[20]*100, color = '#f9c116', marker = '<', alpha = 0.33)
plt.scatter(fezA[21], apl[21]*100, label = 'DVPC -50b Run 1', color = '#f9c116', marker = '>')
plt.scatter(fezB[21], apl[21]*100, color = '#f9c116', marker = '>')
plt.scatter(fezA[22], apl[22]*100, label = 'DVPC -50b Run 2', color = '#f9c116', marker = '>', alpha = 0.66)
plt.scatter(fezB[22], apl[22]*100, color = '#f9c116', marker = '>', alpha = 0.66)
plt.scatter(fezA[23], apl[23]*100, label = 'DVPC -60b Run 1', color = '#f9c116', marker = 'v')
plt.scatter(fezB[23], apl[23]*100, color = '#f9c116', marker = 'v')
plt.scatter(fezA[24], apl[24]*100, label = 'DVPC -60b Run 2', color = '#f9c116', marker = 'v', alpha = 0.66)
plt.scatter(fezB[24], apl[24]*100, color = '#f9c116', marker = 'v', alpha = 0.66)
#plt.scatter(fezA[25], apl[25]*100, label = 'POPC 1b Run 1', color = 'black', marker = '^')
#plt.scatter(fezB[25], apl[25]*100, color = 'black', marker = '^')
#plt.scatter(fezA[26], apl[26]*100, label = 'POPC 1b Run 2', color = 'black', marker = '^', alpha = 0.66)
#plt.scatter(fezB[26], apl[26]*100, color = 'black', marker = '^', alpha = 0.66)
#plt.scatter(fezA[27], apl[27]*100, label = 'POPC -50b Run 1', color = 'black', marker = '<')
#plt.scatter(fezB[27], apl[27]*100, color = 'black', marker = '<')
#plt.scatter(fezA[28], apl[28]*100, label = 'POPC -50b Run 2', color = 'black', marker = '<', alpha = 0.66)
#plt.scatter(fezB[28], apl[28]*100, color = 'black', marker = '<', alpha = 0.66)

v = [0,102,50,120]
plt.axis(v)
plt.axvspan(100, 102, facecolor = 'black', alpha = 0.2)
plt.ylabel("Area Per Lipid ($\AA^2$)", fontsize = '24')
plt.xlabel("Opening Time (ns)", fontsize = '24')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(fontsize = 'small', loc = 'center left', bbox_to_anchor=(1.05,0.5))
plt.title("Area per Lipid vs Time when F < 6.1, E > 11.2, Z > 11.2 ($\AA$)", fontsize = 24)
plt.savefig("aplvsfez-gmm2.png", format='png', dpi=300)
plt.savefig("aplvsfez-gmm2.svg", format='svg', dpi=300)                               
