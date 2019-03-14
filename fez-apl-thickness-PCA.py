import numpy as np
import MDAnalysis
import MDAnalysis.analysis.leaflet
import MDAnalysis.analysis.rms
import matplotlib
import matplotlib.pyplot as plt
import sklearn
from sklearn.decomposition import PCA
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

def discatF (A, B, x, y) :
    U = MDAnalysis.Universe(A, B)
    Fen_A = []
    Fen_B = []

    for ts in U.trajectory:
        Fen_A.append(dis('resid 324 and segid A', 'resid 198 and segid B', U))
        Fen_B.append(dis('resid 324 and segid B', 'resid 198 and segid A', U))
    
    FenA = np.stack((Fen_A, x, y), axis = 1)
    FenB = np.stack((Fen_B, x, y), axis = 1)
    Fen = np.concatenate((FenA, FenB), axis = 0)
    
    return Fen

def discatZ (A, B, x, y) :
    U = MDAnalysis.Universe(A, B)
    Zip_A = []
    Zip_B = []

    for ts in U.trajectory:
        Zip_A.append(dis('resid 326 and segid A', 'resid 237 and segid A', U))
        Zip_B.append(dis('resid 326 and segid B', 'resid 237 and segid B', U))
    
    ZipA = np.stack((Zip_A, x, y), axis = 1)
    ZipB = np.stack((Zip_B, x, y), axis = 1)
    Zip = np.concatenate((ZipA, ZipB), axis = 0)
    
    return Zip

def discatE (A, B, x, y) :
    U = MDAnalysis.Universe(A, B)
    Exp_A = []
    Exp_B = []

    for ts in U.trajectory:
        Exp_A.append(dis('resid 322 and segid A', 'resid 212 and segid A', U))
        Exp_B.append(dis('resid 322 and segid B', 'resid 212 and segid B', U))
    
    ExpA = np.stack((Exp_A, x, y), axis = 1)
    ExpB = np.stack((Exp_B, x, y), axis = 1)
    Exp = np.concatenate((ExpA, ExpB), axis = 0)
    
    return Exp

def draw_vector(v0, v1, ax = None):
    ax = ax or plt.gca()
    arrowprops = dict(arrowstyle='->', linewidth = 2, color = 'black', shrinkA = 0, shrinkB = 0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)

##########
# CALCULATIONS
##########

########## BILAYER THICKNESS

##### DFPC
### 1b

u1 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc')
L1 = MDAnalysis.analysis.leaflet.LeafletFinder(u1, 'name P*')

bilayer1 = []

for ts in u1.trajectory :
	bilayer1.append(bilayer_thickness(u1, L1))

u2 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L2 = MDAnalysis.analysis.leaflet.LeafletFinder(u2, 'name P*')

bilayer2 = []

for ts in u2.trajectory :
	bilayer2.append(bilayer_thickness(u2, L2))

### -30b

u3 = MDAnalysis.Universe('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc')
L3 = MDAnalysis.analysis.leaflet.LeafletFinder(u3, 'name P*')

bilayer3 = []

for ts in u3.trajectory :
	bilayer3.append(bilayer_thickness(u3, L3))

u4 = MDAnalysis.Universe('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc')
L4 = MDAnalysis.analysis.leaflet.LeafletFinder(u4, 'name P*')

bilayer4 = []

for ts in u4.trajectory :
	bilayer4.append(bilayer_thickness(u4, L4))

### -40b

bilayer5 = []

u5 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L5 = MDAnalysis.analysis.leaflet.LeafletFinder(u5, 'name P*') #selects P atoms and puts in two leaflet

for ts in u5.trajectory :
	bilayer5.append(bilayer_thickness(u5, L5))

bilayer6 = []

u6 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L6 = MDAnalysis.analysis.leaflet.LeafletFinder(u6, 'name P*') #selects P atoms and puts in two leaflet

for ts in u6.trajectory :
	bilayer6.append(bilayer_thickness(u6, L6))

bilayer7 = []

u7 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L7 = MDAnalysis.analysis.leaflet.LeafletFinder(u7, 'name P*') #selects P atoms and puts in two leaflet

for ts in u7.trajectory :
	bilayer7.append(bilayer_thickness(u7, L7))

##### DOPC
### 1b

u8 = MDAnalysis.Universe('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc')
L8 = MDAnalysis.analysis.leaflet.LeafletFinder(u8, 'name P*')

bilayer8 = []

for ts in u8.trajectory :
	bilayer8.append(bilayer_thickness(u8, L8))

u9 = MDAnalysis.Universe('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L9 = MDAnalysis.analysis.leaflet.LeafletFinder(u9, 'name P*')

bilayer9 = []

for ts in u9.trajectory :
	bilayer9.append(bilayer_thickness(u9, L9))

### -30b

u10 = MDAnalysis.Universe('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc')
L10 = MDAnalysis.analysis.leaflet.LeafletFinder(u10, 'name P*')

bilayer10 = []

for ts in u10.trajectory :
	bilayer10.append(bilayer_thickness(u10, L10))

u11 = MDAnalysis.Universe('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc')
L11 = MDAnalysis.analysis.leaflet.LeafletFinder(u11, 'name P*')

bilayer11 = []

for ts in u11.trajectory :
	bilayer11.append(bilayer_thickness(u11, L11))

### -40b

bilayer12 = []

u12 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L12 = MDAnalysis.analysis.leaflet.LeafletFinder(u12, 'name P*') #selects P atoms and puts in two leaflet

for ts in u12.trajectory :
	bilayer12.append(bilayer_thickness(u12, L12))

bilayer13 = []

u13 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L13 = MDAnalysis.analysis.leaflet.LeafletFinder(u13, 'name P*') #selects P atoms and puts in two leaflet

for ts in u13.trajectory :
	bilayer13.append(bilayer_thickness(u13, L13))

bilayer14 = []

u14 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L14 = MDAnalysis.analysis.leaflet.LeafletFinder(u14, 'name P*') #selects P atoms and puts in two leaflet

for ts in u14.trajectory :
	bilayer14.append(bilayer_thickness(u14, L14))

### -50b

u15 = MDAnalysis.Universe('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L15 = MDAnalysis.analysis.leaflet.LeafletFinder(u15, 'name P*')

bilayer15 = []

for ts in u15.trajectory :
	bilayer15.append(bilayer_thickness(u15, L15))

u16 = MDAnalysis.Universe('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L16 = MDAnalysis.analysis.leaflet.LeafletFinder(u16, 'name P*')

bilayer16 = []

for ts in u16.trajectory :
	bilayer16.append(bilayer_thickness(u16, L16))

##### DVPC
### 1b

u17 = MDAnalysis.Universe('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc')
L17 = MDAnalysis.analysis.leaflet.LeafletFinder(u17, 'name P*')

bilayer17 = []

for ts in u17.trajectory :
	bilayer17.append(bilayer_thickness(u17, L17))

u18 = MDAnalysis.Universe('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
L18 = MDAnalysis.analysis.leaflet.LeafletFinder(u18, 'name P*')

bilayer18 = []

for ts in u18.trajectory :
	bilayer18.append(bilayer_thickness(u18, L18))

### -40b

bilayer19 = []

u19 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
L19 = MDAnalysis.analysis.leaflet.LeafletFinder(u19, 'name P*') #selects P atoms and puts in two leaflet

for ts in u19.trajectory :
	bilayer19.append(bilayer_thickness(u19, L19))

bilayer20 = []

u20 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
L20 = MDAnalysis.analysis.leaflet.LeafletFinder(u20, 'name P*') #selects P atoms and puts in two leaflet

for ts in u20.trajectory :
	bilayer20.append(bilayer_thickness(u20, L20))

bilayer21 = []

u21 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc')
L21 = MDAnalysis.analysis.leaflet.LeafletFinder(u21, 'name P*') #selects P atoms and puts in two leaflet

for ts in u21.trajectory :
	bilayer21.append(bilayer_thickness(u21, L21))

### -50b

u22 = MDAnalysis.Universe('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L22 = MDAnalysis.analysis.leaflet.LeafletFinder(u22, 'name P*')

bilayer22 = []

for ts in u22.trajectory :
	bilayer22.append(bilayer_thickness(u22, L22))

u23 = MDAnalysis.Universe('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L23 = MDAnalysis.analysis.leaflet.LeafletFinder(u23, 'name P*')

bilayer23 = []

for ts in u23.trajectory :
	bilayer23.append(bilayer_thickness(u23, L23))

### -60b

u24 = MDAnalysis.Universe('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc')
L24 = MDAnalysis.analysis.leaflet.LeafletFinder(u24, 'name P*')

bilayer24 = []

for ts in u24.trajectory :
	bilayer24.append(bilayer_thickness(u24, L24))

u25 = MDAnalysis.Universe('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc')
L25 = MDAnalysis.analysis.leaflet.LeafletFinder(u25, 'name P*')

bilayer25 = []

for ts in u25.trajectory :
	bilayer25.append(bilayer_thickness(u25, L25))

##### POPC
### 1b

u26 = MDAnalysis.Universe('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc')
L26 = MDAnalysis.analysis.leaflet.LeafletFinder(u26, 'name P*')

bilayer26 = []

for ts in u26.trajectory :
	bilayer26.append(bilayer_thickness(u26, L26))

u27 = MDAnalysis.Universe('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc')
L27 = MDAnalysis.analysis.leaflet.LeafletFinder(u27, 'name P*')

bilayer27 = []

for ts in u27.trajectory :
	bilayer27.append(bilayer_thickness(u27, L27))

### -50b

u28 = MDAnalysis.Universe('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
L28 = MDAnalysis.analysis.leaflet.LeafletFinder(u28, 'name P*')

bilayer28 = []

for ts in u28.trajectory :
	bilayer28.append(bilayer_thickness(u28, L28))

u29 = MDAnalysis.Universe('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc')
L29 = MDAnalysis.analysis.leaflet.LeafletFinder(u29, 'name P*')

bilayer29 = []

for ts in u29.trajectory :
	bilayer29.append(bilayer_thickness(u29, L29))

########## APL

##### DFPC
### 1b

d1 = np.loadtxt('../../DFPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)

d2 = np.loadtxt('../../DFPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)

### -30b

d3 = np.loadtxt('../../DFPC/Backwards/md/30bar/repeat1/analysis/apl.xvg', skiprows = 15)

d4 = np.loadtxt('../../DFPC/Backwards/md/30bar/repeat2/analysis/apl.xvg', skiprows = 15)

### -40b

d5 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)

d6 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)

d7 = np.loadtxt('../../DFPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)

##### DOPC
### 1b

d8 = np.loadtxt('../../DOPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)

d9 = np.loadtxt('../../DOPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)

### -30b

d10 = np.loadtxt('../../DOPC/Backwards/md/30bar/repeat1/analysis/apl.xvg', skiprows = 15)

d11 = np.loadtxt('../../DOPC/Backwards/md/30bar/repeat2/analysis/apl.xvg', skiprows = 15)

### -40b

d12 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)

d13 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)

d14 = np.loadtxt('../../DOPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)

### -50b

d15 = np.loadtxt('../../DOPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)

d16 = np.loadtxt('../../DOPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)

##### DVPC
### 1b

d17 = np.loadtxt('../../DVPC/Backwards/md/1bar/analysis/apl.xvg', skiprows = 15)

d18 = np.loadtxt('../../DVPC/Backwards/md/1bar/repeat/analysis/apl.xvg', skiprows = 15)

### -40b

d19 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat1/analysis/apl.xvg', skiprows = 15)

d20 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat2/analysis/apl.xvg', skiprows = 15)

d21 = np.loadtxt('../../DVPC/Backwards/md/40bar/repeat3/analysis/apl.xvg', skiprows = 15)

### -50b

d22 = np.loadtxt('../../DVPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)

d23 = np.loadtxt('../../DVPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)

### -60b

d24 = np.loadtxt('../../DVPC/Backwards/md/60bar/repeat1/analysis/apl.xvg', skiprows = 15)

d25 = np.loadtxt('../../DVPC/Backwards/md/60bar/repeat2/analysis/apl.xvg', skiprows = 15)

##### POPC
### 1b

d26 = np.loadtxt('../../POPC/Backwards/md/1bar/repeat1/analysis/apl.xvg', skiprows = 15)

d27 = np.loadtxt('../../POPC/Backwards/md/1bar/repeat2/analysis/apl.xvg', skiprows = 15)

### -50b

d28 = np.loadtxt('../../POPC/Backwards/md/50bar/repeat1/analysis/apl.xvg', skiprows = 15)

d29 = np.loadtxt('../../POPC/Backwards/md/50bar/repeat2/analysis/apl.xvg', skiprows = 15)

########## FENESTRATION

FenTot = np.zeros((1,3))

##### DFPC
### 1b

FenTot = np.concatenate((FenTot, discatF('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer1, d1[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer2, d2[:,1])), axis=0)

### -30b

FenTot = np.concatenate((FenTot, discatF('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', bilayer3, d3[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', bilayer4, d4[:,1])), axis=0)

### -40b
FenTot = np.concatenate((FenTot, discatF('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer5, d5[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer6, d6[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer7, d7[:,1])), axis=0)

##### DOPC
### 1b

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer8, d8[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer9, d9[:,1])), axis=0)

### -30b

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', bilayer10, d10[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', bilayer11, d11[:,1])), axis=0)

### -40b

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer12, d12[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer13, d13[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer14, d14[:,1])), axis=0)

### -50b

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer15, d15[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer16, d16[:,1])), axis=0)

##### DVPC
### 1b
FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer17, d17[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer18, d18[:,1])), axis=0)

### -40b

FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer19, d19[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer20, d20[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer21, d21[:,1])), axis=0)

### -50b

FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer22, d22[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer23, d23[:,1])), axis=0)

### -60b

FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', bilayer24, d24[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', bilayer25, d25[:,1])), axis=0)

##### POPC
### 1b
FenTot = np.concatenate((FenTot, discatF('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc', bilayer26, d26[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc', bilayer27, d27[:,1])), axis=0)

### -50b

FenTot = np.concatenate((FenTot, discatF('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer28, d28[:,1])), axis=0)

FenTot = np.concatenate((FenTot, discatF('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer29, d29[:,1])), axis=0)

########## ZIPPER

ZipTot = np.zeros((1,3))

##### DFPC
### 1b

ZipTot = np.concatenate((ZipTot, discatZ('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer1, d1[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer2, d2[:,1])), axis=0)

### -30b

ZipTot = np.concatenate((ZipTot, discatZ('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', bilayer3, d3[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', bilayer4, d4[:,1])), axis=0)

### -40b
ZipTot = np.concatenate((ZipTot, discatZ('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer5, d5[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer6, d6[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer7, d7[:,1])), axis=0)

##### DOPC
### 1b

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer8, d8[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer9, d9[:,1])), axis=0)

### -30b

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', bilayer10, d10[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', bilayer11, d11[:,1])), axis=0)

### -40b

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer12, d12[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer13, d13[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer14, d14[:,1])), axis=0)

### -50b

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer15, d15[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer16, d16[:,1])), axis=0)

##### DVPC
### 1b
ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer17, d17[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer18, d18[:,1])), axis=0)

### -40b

ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer19, d19[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer20, d20[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer21, d21[:,1])), axis=0)

### -50b

ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer22, d22[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer23, d23[:,1])), axis=0)

### -60b

ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', bilayer24, d24[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', bilayer25, d25[:,1])), axis=0)

##### POPC
### 1b
ZipTot = np.concatenate((ZipTot, discatZ('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc', bilayer26, d26[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc', bilayer27, d27[:,1])), axis=0)

### -50b

ZipTot = np.concatenate((ZipTot, discatZ('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer28, d28[:,1])), axis=0)

ZipTot = np.concatenate((ZipTot, discatZ('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer29, d29[:,1])), axis=0)

########## EXPANSION

ExpTot = np.zeros((1,3))

##### DFPC
### 1b

ExpTot = np.concatenate((ExpTot, discatE('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer1, d1[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer2, d2[:,1])), axis=0)

### -30b

ExpTot = np.concatenate((ExpTot, discatE('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', bilayer3, d3[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', bilayer4, d4[:,1])), axis=0)

### -40b
ExpTot = np.concatenate((ExpTot, discatE('../../DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer5, d5[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer6, d6[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer7, d7[:,1])), axis=0)

##### DOPC
### 1b

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer8, d8[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer9, d9[:,1])), axis=0)

### -30b

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', bilayer10, d10[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', bilayer11, d11[:,1])), axis=0)

### -40b

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer12, d12[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer13, d13[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer14, d14[:,1])), axis=0)

### -50b

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer15, d15[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer16, d16[:,1])), axis=0)

##### DVPC
### 1b
ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc', bilayer17, d17[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', bilayer18, d18[:,1])), axis=0)

### -40b

ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', bilayer19, d19[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', bilayer20, d20[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', bilayer21, d21[:,1])), axis=0)

### -50b

ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer22, d22[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer23, d23[:,1])), axis=0)

### -60b

ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', bilayer24, d24[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', bilayer25, d25[:,1])), axis=0)

##### POPC
### 1b
ExpTot = np.concatenate((ExpTot, discatE('../../POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc', bilayer26, d26[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc', bilayer27, d27[:,1])), axis=0)

### -50b

ExpTot = np.concatenate((ExpTot, discatE('../../POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', bilayer28, d28[:,1])), axis=0)

ExpTot = np.concatenate((ExpTot, discatE('../../POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '../../POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', bilayer29, d29[:,1])), axis=0)

########## PCA

X = FenTot[1:, [1,0]]
pca1a = PCA(n_components=2)
pca1a.fit(X)
pca1b = PCA(n_components=1)
pca1b.fit(X)
X_trans = pca1b.transform(X)
X_new = pca1b.inverse_transform(X_trans)

Y = ZipTot[1:, [1,0]]
pca2a = PCA(n_components=2)
pca2a.fit(Y)
pca2b = PCA(n_components=1)
pca2b.fit(Y)
Y_trans = pca2b.transform(Y)
Y_new = pca2b.inverse_transform(Y_trans)

Z = ExpTot[1:, [1,0]]
pca3a = PCA(n_components=2)
pca3a.fit(Z)
pca3b = PCA(n_components=1)
pca3b.fit(Z)
Z_trans = pca3b.transform(Z)
Z_new = pca3b.inverse_transform(Z_trans)

x = FenTot[1:, [2,0]]
pca4a = PCA(n_components=2)
pca4a.fit(x)
pca4b = PCA(n_components=1)
pca4b.fit(x)

x_trans = pca4b.transform(x)
x_new = pca4b.inverse_transform(x_trans)

y = ZipTot[1:, [2,0]]
pca5a = PCA(n_components=2)
pca5a.fit(y)
pca5b = PCA(n_components=1)
pca5b.fit(y)
y_trans = pca5b.transform(y)
y_new = pca5b.inverse_transform(y_trans)

z = ExpTot[1:, [2,0]]
pca6a = PCA(n_components=2)
pca6a.fit(z)
pca6b = PCA(n_components=1)
pca6b.fit(z)
z_trans = pca6b.transform(z)
z_new = pca6b.inverse_transform(z_trans)

##########
# PLOTTING
##########

plt.figure(figsize = (18,6))

plt.subplot(131)
plt.scatter(X[0:1414,0], X[0:1414:,1], alpha = 0.4, marker = 'o', color = 'black')
plt.scatter(X[1414:3232,0], X[1414:3232,1], alpha = 0.4, marker = 'o', color = 'grey')
plt.scatter(X[3232:5050,0], X[3232:5050,1], alpha = 0.4, marker = 'o', color = 'blue')
for length, vector in zip(pca1a.explained_variance_, pca1a.components_):
    v = vector * 3 * np.sqrt(length)
    draw_vector(pca1a.mean_, pca1a.mean_ + v)
#plt.scatter(X_new[:,0], X_new[:,1], alpha = 0.8, marker = 'o', color = 'blue')
plt.ylabel('Distance ($\AA$)', fontsize = 14)
plt.xlabel('Bilayer Thickness ($\AA$)', fontsize = 14)
q = [28, 44, 4, 18]
plt.axis(q)
plt.title('Fenestration', fontsize = 18)

plt.subplot(132)
plt.scatter(Y[0:1414,0], Y[0:1414,1], alpha = 0.4, color = 'black', marker = '^')
plt.scatter(Y[1414:3232,0], Y[1414:3232,1], alpha = 0.4, color = 'grey', marker = '^')
plt.scatter(Y[3232:5050,0], Y[3232:5050,1], alpha = 0.4, color = 'green', marker = '^')
for length, vector in zip(pca2a.explained_variance_, pca2a.components_):
    v = vector * 3 * np.sqrt(length)
    draw_vector(pca2a.mean_, pca2a.mean_ + v)
#plt.scatter(Y_new[:,0], Y_new[:,1], alpha = 0.8, color = 'green', marker = '^')
plt.ylabel('Distance ($\AA$)', fontsize = 14)
plt.xlabel('Bilayer Thickness ($\AA$)', fontsize = 14)
plt.axis(q)
plt.title('Zipper', fontsize = 18)
    
plt.subplot(133)
plt.scatter(Z[0:1414,0], Z[0:1414,1], alpha = 0.4, marker = 'D', color = 'black')
plt.scatter(Z[1414:3232,0], Z[1414:3232,1], alpha = 0.4, marker = 'D', color = 'grey') 
plt.scatter(Z[3232:5050,0], Z[3232:5050,1], alpha = 0.4, marker = 'D', color = 'red') 
for length, vector in zip(pca3a.explained_variance_, pca3a.components_):
    v = vector * 3 * np.sqrt(length)
    draw_vector(pca3a.mean_, pca3a.mean_ + v)
#plt.scatter(Z_new[:,0], Z_new[:,1], alpha = 0.8, marker = 'D', color = 'red')
plt.ylabel('Distance ($\AA$)', fontsize = 14)
plt.xlabel('Bilayer Thickness ($\AA$)', fontsize = 14)
plt.axis(q)
plt.title('Expansion', fontsize = 18)

plt.savefig("alldist-thickness.png", format='png', dpi=300)
plt.savefig("alldist-thickness.svg", format='svg', dpi=300)  

plt.clf()

plt.figure(figsize = (18,6))

plt.subplot(131)
plt.scatter(x[0:1414,0], x[0:1414:,1], alpha = 0.4, marker = 'o', color = 'black')
plt.scatter(x[1414:3232,0], x[1414:3232,1], alpha = 0.4, marker = 'o', color = 'grey')
plt.scatter(x[3232:5050,0], x[3232:5050,1], alpha = 0.4, marker = 'o', color = 'blue')
for length, vector in zip(pca4a.explained_variance_, pca4a.components_):
    v = vector * 3 * np.sqrt(length)
    draw_vector(pca4a.mean_, pca4a.mean_ + v)
#plt.scatter(x_new[:,0], x_new[:,1], alpha = 0.8, marker = 'o', color = 'blue')
plt.ylabel('Distance ($\AA$)', fontsize = 14)
plt.xlabel('Area Per Lipid ($\AA^2$)', fontsize = 14)
w = [0.6, 1.2, 4, 18]
plt.axis(w)
plt.title('Fenestration', fontsize = 18)

plt.subplot(132)
plt.scatter(y[0:1414,0], y[0:1414,1], alpha = 0.4, color = 'black', marker = '^')
plt.scatter(y[1414:3232,0], y[1414:3232,1], alpha = 0.4, color = 'grey', marker = '^')
plt.scatter(y[3232:5050,0], y[3232:5050,1], alpha = 0.4, color = 'green', marker = '^')
for length, vector in zip(pca5a.explained_variance_, pca5a.components_):
    v = vector * 3 * np.sqrt(length)
    draw_vector(pca5a.mean_, pca5a.mean_ + v)
#plt.scatter(y_new[:,0], y_new[:,1], alpha = 0.8, color = 'green', marker = '^')
plt.ylabel('Distance ($\AA$)', fontsize = 14)
plt.xlabel('Area Per Lipid ($\AA^2$)', fontsize = 14)
plt.axis(w)
plt.title('Zipper', fontsize = 18)
    
plt.subplot(133)
plt.scatter(z[0:1414,0], z[0:1414,1], alpha = 0.4, marker = 'D', color = 'black')
plt.scatter(z[1414:3232,0], z[1414:3232,1], alpha = 0.4, marker = 'D', color = 'grey') 
plt.scatter(z[3232:5050,0], z[3232:5050,1], alpha = 0.4, marker = 'D', color = 'red')  
for length, vector in zip(pca6a.explained_variance_, pca6a.components_):
    v = vector * 3 * np.sqrt(length)
    draw_vector(pca6a.mean_, pca6a.mean_ + v)
#plt.scatter(z_new[:,0], z_new[:,1], alpha = 0.8, marker = 'D', color = 'red')
plt.ylabel('Distance ($\AA$)', fontsize = 14)
plt.xlabel('Area Per Lipid ($\AA^2$)', fontsize = 14)
plt.axis(w)
plt.title('Expansion', fontsize = 18)

plt.savefig("alldist-apl.png", format='png', dpi=300)
plt.savefig("alldist-apl.svg", format='svg', dpi=300)
