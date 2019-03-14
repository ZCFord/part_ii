import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import sklearn
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture as GMM
import numpy as np
from matplotlib.patches import Ellipse

rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18

##########
# DEFINITIONS
##########

def dis (sel1, sel2, x) :
	res_i = x.select_atoms(sel1 + ' and name CA').positions
	res_j = x.select_atoms(sel2 + ' and name CA').positions
	dis = np.linalg.norm(res_i - res_j)
	return dis

def F_Z(A, B, C) :
        u = MDAnalysis.Universe(A, B)
        protein = u.select_atoms('protein')
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
        FZ_a = np.stack((Zip_A, Fen_A), axis = 1)
        FZ_b = np.stack((Zip_B, Fen_B), axis = 1)
        C.append(np.concatenate((FZ_a, FZ_b), axis = 0))
        C = np.squeeze(C)
        return C

def F_E(A, B, C) :
        u = MDAnalysis.Universe(A, B)
        protein = u.select_atoms('protein')
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
        FE_a = np.stack((Exp_A, Fen_A), axis = 1)
        FE_b = np.stack((Exp_B, Fen_B), axis = 1)
        C.append(np.concatenate((FE_a, FE_b), axis = 0))
        C = np.squeeze(C)
        return C

def draw_ellipse(position, covariance, ax=None, **kwargs):
    ax = ax or plt.gca()
    
    if covariance.shape == (2,2):
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1,0], U[0,0]))
        width, height = 2 * np.sqrt(s)
    else :
        angle = 0
        width, height = 2 * np.sqrt(covariance)
        
    for nsig in range(1,4):
        ax.add_patch(Ellipse(position, nsig * width, nsig * height, angle, **kwargs))
        
def plot_gmm(gmm, X, label = True, ax = None):
    w_factor = 0.2 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
        draw_ellipse(pos, covar, alpha = w * w_factor, edgecolor = 'black', facecolor = 'yellow')

def counting(A):
    count = 0
    for i in A[:,3]:
        if i > 0.5:
            count = count + 1
    percent = (float(count)/202)*100
    return percent

def counting_harsh(A):
    count = 0
    for i in A[:,3]:
        if i > 0.9:
            count = count + 1
    percent = (float(count)/202)*100
    return percent

##########
# STRUCTURES
##########

up = MDAnalysis.Universe('/sansom/s150/pemb4066/Documents/PART_II/USEFUL/structures/UP_correctpro.pdb')
up_pro = up.select_atoms('protein')

dn = MDAnalysis.Universe('/sansom/s150/pemb4066/Documents/PART_II/USEFUL/structures/DN_correctResid.pdb')
dn_pro = dn.select_atoms('protein')

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
# CALCULATIONS
##########

####### FZ

##### DFPC
### 1b

F1b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/analysis/md-c.xtc', F1b1FZ)
F1b1FZ = np.squeeze(F1b1FZ)

F1b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', F1b2FZ)
F1b2FZ = np.squeeze(F1b2FZ)

### -30b
F30b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', F30b1FZ)
F30b1FZ = np.squeeze(F30b1FZ)

F30b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', F30b2FZ)
F30b2FZ = np.squeeze(F30b2FZ)

### -40b

F40b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', F40b1FZ)
F40b1FZ = np.squeeze(F40b1FZ)

F40b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', F40b2FZ)
F40b2FZ = np.squeeze(F40b2FZ)

F40b3FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', F40b3FZ)
F40b3FZ = np.squeeze(F40b3FZ)

##### DOPC
### 1b

O1b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/1bar/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/1bar/analysis/md-c.xtc', O1b1FZ)
O1b1FZ = np.squeeze(O1b1FZ)

O1b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', O1b2FZ)
O1b2FZ = np.squeeze(O1b2FZ)

### -30b

O30b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', O30b1FZ)
O30b1FZ = np.squeeze(O30b1FZ)

O30b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', O30b2FZ)
O30b2FZ = np.squeeze(O30b2FZ)

### -40b

O40b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', O40b1FZ)
O40b1FZ = np.squeeze(O40b1FZ)

O40b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', O40b2FZ)
O40b2FZ = np.squeeze(O40b2FZ)

O40b3FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', O40b3FZ)
O40b3FZ = np.squeeze(O40b3FZ)

### -50b

O50b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', O50b1FZ)
O50b1FZ = np.squeeze(O50b1FZ)

O50b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', O50b2FZ)
O50b2FZ = np.squeeze(O50b2FZ)

##### DVPC
### 1b

V1b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/1bar/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/1bar/analysis/md-c.xtc', V1b1FZ)
V1b1FZ = np.squeeze(V1b1FZ)

V1b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', V1b2FZ)
V1b2FZ = np.squeeze(V1b2FZ)

### -40b

V40b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', V40b1FZ)
V40b1FZ = np.squeeze(V40b1FZ)

V40b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', V40b2FZ)
V40b2FZ = np.squeeze(V40b2FZ)

V40b3FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', V40b3FZ)
V40b3FZ = np.squeeze(V40b3FZ)

### -50b

V50b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', V50b1FZ)
V50b1FZ = np.squeeze(V50b1FZ)

V50b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', V50b2FZ)
V50b2FZ = np.squeeze(V50b2FZ)

### -60b

V60b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', V60b1FZ)
V60b1FZ = np.squeeze(V60b1FZ)

V60b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', V60b2FZ)
V60b2FZ = np.squeeze(V60b2FZ)

####### FE

##### DFPC
### 1b

F1b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/analysis/md-c.xtc', F1b1FE)
F1b1FE = np.squeeze(F1b1FE)

F1b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', F1b2FE)
F1b2FE = np.squeeze(F1b2FE)

### -30b
F30b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', F30b1FE)
F30b1FE = np.squeeze(F30b1FE)

F30b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', F30b2FE)
F30b2FE = np.squeeze(F30b2FE)

### -40b

F40b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', F40b1FE)
F40b1FE = np.squeeze(F40b1FE)

F40b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', F40b2FE)
F40b2FE = np.squeeze(F40b2FE)

F40b3FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', F40b3FE)
F40b3FE = np.squeeze(F40b3FE)

##### DOPC
### 1b

O1b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/1bar/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/1bar/analysis/md-c.xtc', O1b1FE)
O1b1FE = np.squeeze(O1b1FE)

O1b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', O1b2FE)
O1b2FE = np.squeeze(O1b2FE)

### -30b

O30b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc', O30b1FE)
O30b1FE = np.squeeze(O30b1FE)

O30b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/30bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/30bar/repeat2/analysis/md-c.xtc', O30b2FE)
O30b2FE = np.squeeze(O30b2FE)

### -40b

O40b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', O40b1FE)
O40b1FE = np.squeeze(O40b1FE)

O40b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', O40b2FE)
O40b2FE = np.squeeze(O40b2FE)

O40b3FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', O40b3FE)
O40b3FE = np.squeeze(O40b3FE)

### -50b

O50b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', O50b1FE)
O50b1FE = np.squeeze(O50b1FE)

O50b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DOPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', O50b2FE)
O50b2FE = np.squeeze(O50b2FE)

##### DVPC
### 1b

V1b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/1bar/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/1bar/analysis/md-c.xtc', V1b1FE)
V1b1FE = np.squeeze(V1b1FE)

V1b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc', V1b2FE)
V1b2FE = np.squeeze(V1b2FE)

### -40b

V40b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc', V40b1FE)
V40b1FE = np.squeeze(V40b1FE)

V40b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc', V40b2FE)
V40b2FE = np.squeeze(V40b2FE)

V40b3FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat3/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/40bar/repeat3/analysis/md-c.xtc', V40b3FE)
V40b3FE = np.squeeze(V40b3FE)

### -50b

V50b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', V50b1FE)
V50b1FE = np.squeeze(V50b1FE)

V50b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', V50b2FE)
V50b2FE = np.squeeze(V50b2FE)

### -60b

V60b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc', V60b1FE)
V60b1FE = np.squeeze(V60b1FE)

V60b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/60bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DVPC/Backwards/md/60bar/repeat2/analysis/md-c.xtc', V60b2FE)
V60b2FE = np.squeeze(V60b2FE)

##########
# CLUSTERING
##########

FZ_tot = np.concatenate((F1b1FZ, F1b2FZ, F30b1FZ, F30b2FZ, F40b1FZ, F40b2FZ, F40b3FZ, O1b1FZ, O1b2FZ, O30b1FZ, O30b2FZ, O40b1FZ, O40b2FZ, O40b3FZ, O50b1FZ, O50b2FZ, V1b1FZ, V1b2FZ, V40b1FZ, V40b2FZ, V40b3FZ, V50b1FZ, V50b2FZ, V60b1FZ, V60b2FZ), axis = 0)

k_means = KMeans(n_clusters = 3)
k_means.fit(FZ_tot)
y_k_means = k_means.predict(FZ_tot)

gmm = GMM(n_components = 2, covariance_type='full')
labels = gmm.fit(FZ_tot).predict(FZ_tot)

FE_tot = np.concatenate((F1b1FE, F1b2FE, F30b1FE, F30b2FE, F40b1FE, F40b2FE, F40b3FE, O1b1FE, O1b2FE, O30b1FE, O30b2FE, O40b1FE, O40b2FE, O40b3FE, O50b1FE, O50b2FE, V1b1FE, V1b2FE, V40b1FE, V40b2FE, V40b3FE, V50b1FE, V50b2FE, V60b1FE, V60b2FE), axis = 0)

kmeans = KMeans(n_clusters = 3)
kmeans.fit(FE_tot)
y_kmeans = kmeans.predict(FE_tot)

Gmm = GMM(n_components = 2, covariance_type='full')
label = Gmm.fit(FE_tot).predict(FE_tot)

##########
# BAR CHART
##########

####### FZ

probs = gmm.predict_proba(FZ_tot)

probFZ = np.concatenate((FZ_tot, probs), axis = 1)

F1b1z = probFZ[0:202,:]
F1b2z = probFZ[202:404,:]
F30b1z = probFZ[404:606,:]
F30b2z = probFZ[606:808,:]
F40b1z = probFZ[808:1010,:]
F40b2z = probFZ[1010:1212,:]
F40b3z = probFZ[1212:1414,:]
O1b1z = probFZ[1414:1616,:]
O1b2z = probFZ[1616:1818,:]
O30b1z = probFZ[1818:2020,:]
O30b2z = probFZ[2020:2222,:]
O40b1z = probFZ[2222:2424,:]
O40b2z = probFZ[2424:2626,:]
O40b3z = probFZ[2626:2828,:]
O50b1z = probFZ[2828:3030,:]
O50b2z = probFZ[3030:3232,:]
V1b1z = probFZ[3232:3434,:]
V1b2z = probFZ[3434:3636,:]
V40b1z = probFZ[3636:3838,:]
V40b2z = probFZ[3838:4040,:]
V40b3z = probFZ[4040:4242,:]
V50b1z = probFZ[4242:4444,:]
V50b2z = probFZ[4444:4646,:]
V60b1z = probFZ[4646:4848,:]
V60b2z = probFZ[4848:5050,:]

zpercentages_harsh = ((counting_harsh(F1b1z)), (counting_harsh(F1b2z)), (counting_harsh(F30b1z)), (counting_harsh(F30b2z)), (counting_harsh(F40b1z)), (counting_harsh(F40b2z)), (counting_harsh(F40b3z)), (counting_harsh(O1b1z)), (counting_harsh(O1b2z)), (counting_harsh(O30b1z)), (counting_harsh(O30b2z)), (counting_harsh(O40b1z)), (counting_harsh(O40b2z)), (counting_harsh(O40b3z)), (counting_harsh(O50b1z)), (counting_harsh(O50b2z)), (counting_harsh(V1b1z)), (counting_harsh(V1b2z)), (counting_harsh(V40b1z)), (counting_harsh(V40b2z)), (counting_harsh(V40b3z)), (counting_harsh(V50b1z)), (counting_harsh(V50b2z)), (counting_harsh(V60b1z)), (counting_harsh(V60b2z)))

zpercentages_diff = ((counting(F1b1z) - counting_harsh(F1b1z)), (counting(F1b2z) - counting_harsh(F1b2z)), (counting(F30b1z) - counting_harsh(F30b1z)), (counting(F30b2z) - counting_harsh(F30b2z)), (counting(F40b1z) - counting_harsh(F40b1z)), (counting(F40b2z) - counting_harsh(F40b2z)), (counting(F40b3z) - counting_harsh(F40b3z)), (counting(O1b1z) - counting_harsh(O1b1z)), (counting(O1b2z) - counting_harsh(O1b2z)), (counting(O30b1z) - counting_harsh(O30b1z)), (counting(O30b2z) - counting_harsh(O30b2z)), (counting(O40b1z) - counting_harsh(O40b1z)), (counting(O40b2z) - counting_harsh(O40b2z)), (counting(O40b3z) - counting_harsh(O40b3z)), (counting(O50b1z) - counting_harsh(O50b1z)), (counting(O50b2z) - counting_harsh(O50b2z)), (counting(V1b1z) - counting_harsh(V1b1z)), (counting(V1b2z) - counting_harsh(V1b2z)), (counting(V40b1z) - counting_harsh(V40b1z)), (counting(V40b2z) - counting_harsh(V40b2z)), (counting(V40b3z) - counting_harsh(V40b3z)), (counting(V50b1z) - counting_harsh(V50b1z)), (counting(V50b2z) - counting_harsh(V50b2z)), (counting(V60b1z) - counting_harsh(V60b1z)), (counting(V60b2z) - counting_harsh(V60b2z)))

F1bz = (counting(F1b1z)+counting(F1b2z))/2
F30bz = (counting(F30b1z)+counting(F30b2z))/2
F40bz = (counting(F40b1z)+counting(F40b2z)+counting(F40b3z))/3
O1bz = (counting(O1b1z)+counting(O1b2z))/2
O30bz = (counting(O30b1z)+counting(O30b2z))/2
O40bz = (counting(O40b1z)+counting(O40b2z)+counting(O40b3z))/3
O50bz = (counting(O50b1z)+counting(O50b2z))/2
V1bz = (counting(V1b1z)+counting(V1b2z))/2
V40bz = (counting(V40b1z)+counting(V40b2z)+counting(V40b3z))/3
V50bz = (counting(V50b1z)+counting(V50b2z))/2
V60bz = (counting(V60b1z)+counting(V60b2z))/2

F1bhz = (counting_harsh(F1b1z)+counting_harsh(F1b2z))/2
F30bhz = (counting_harsh(F30b1z)+counting_harsh(F30b2z))/2
F40bhz = (counting_harsh(F40b1z)+counting_harsh(F40b2z)+counting_harsh(F40b3z))/3
O1bhz = (counting_harsh(O1b1z)+counting_harsh(O1b2z))/2
O30bhz = (counting_harsh(O30b1z)+counting_harsh(O30b2z))/2
O40bhz = (counting_harsh(O40b1z)+counting_harsh(O40b2z)+counting_harsh(O40b3z))/3
O50bhz = (counting_harsh(O50b1z)+counting_harsh(O50b2z))/2
V1bhz = (counting_harsh(V1b1z)+counting_harsh(V1b2z))/2
V40bhz = (counting_harsh(V40b1z)+counting_harsh(V40b2z)+counting_harsh(V40b3z))/3
V50bhz = (counting_harsh(V50b1z)+counting_harsh(V50b2z))/2
V60bhz = (counting_harsh(V60b1z)+counting_harsh(V60b2z))/2

F1bdz = F1bz - F1bhz
F30bdz = F30bz - F30bhz
F40bdz = F40bz - F40bhz
O1bdz = O1bz - O1bhz
O30bdz = O30bz - O30bhz
O40bdz = O40bz - O40bhz
O50bdz = O50bz - O50bhz
V1bdz = V1bz - V1bhz
V40bdz = V40bz - V40bhz
V50bdz = V50bz - V50bhz
V60bdz = V60bz - V60bhz

zpercents_harsh = (F1bhz, F30bhz, F40bhz, O1bhz, O30bhz, O40bhz, O50bhz, V1bhz, V40bhz, V50bhz, V60bhz)
zpercents_dff = (F1bdz, F30bdz, F40bdz, O1bdz, O30bdz, O40bdz, O50bdz, V1bdz, V40bdz, V50bdz, V60bdz)

####### FE

prob = Gmm.predict_proba(FE_tot)

probFE = np.concatenate((FE_tot, prob), axis = 1)

F1b1e = probFE[0:202,:]
F1b2e = probFE[202:404,:]
F30b1e = probFE[404:606,:]
F30b2e = probFE[606:808,:]
F40b1e = probFE[808:1010,:]
F40b2e = probFE[1010:1212,:]
F40b3e = probFE[1212:1414,:]
O1b1e = probFE[1414:1616,:]
O1b2e = probFE[1616:1818,:]
O30b1e = probFE[1818:2020,:]
O30b2e = probFE[2020:2222,:]
O40b1e = probFE[2222:2424,:]
O40b2e = probFE[2424:2626,:]
O40b3e = probFE[2626:2828,:]
O50b1e = probFE[2828:3030,:]
O50b2e = probFE[3030:3232,:]
V1b1e = probFE[3232:3434,:]
V1b2e = probFE[3434:3636,:]
V40b1e = probFE[3636:3838,:]
V40b2e = probFE[3838:4040,:]
V40b3e = probFE[4040:4242,:]
V50b1e = probFE[4242:4444,:]
V50b2e = probFE[4444:4646,:]
V60b1e = probFE[4646:4848,:]
V60b2e = probFE[4848:5050,:]

epercentages_harsh = ((counting_harsh(F1b1e)), (counting_harsh(F1b2e)), (counting_harsh(F30b1e)), (counting_harsh(F30b2e)), (counting_harsh(F40b1e)), (counting_harsh(F40b2e)), (counting_harsh(F40b3e)), (counting_harsh(O1b1e)), (counting_harsh(O1b2e)), (counting_harsh(O30b1e)), (counting_harsh(O30b2e)), (counting_harsh(O40b1e)), (counting_harsh(O40b2e)), (counting_harsh(O40b3e)), (counting_harsh(O50b1e)), (counting_harsh(O50b2e)), (counting_harsh(V1b1e)), (counting_harsh(V1b2e)), (counting_harsh(V40b1e)), (counting_harsh(V40b2e)), (counting_harsh(V40b3e)), (counting_harsh(V50b1e)), (counting_harsh(V50b2e)), (counting_harsh(V60b1e)), (counting_harsh(V60b2e)))

epercentages_diff = ((counting(F1b1e) - counting_harsh(F1b1e)), (counting(F1b2e) - counting_harsh(F1b2e)), (counting(F30b1e) - counting_harsh(F30b1e)), (counting(F30b2e) - counting_harsh(F30b2e)), (counting(F40b1e) - counting_harsh(F40b1e)), (counting(F40b2e) - counting_harsh(F40b2e)), (counting(F40b3e) - counting_harsh(F40b3e)), (counting(O1b1e) - counting_harsh(O1b1e)), (counting(O1b2e) - counting_harsh(O1b2e)), (counting(O30b1e) - counting_harsh(O30b1e)), (counting(O30b2e) - counting_harsh(O30b2e)), (counting(O40b1e) - counting_harsh(O40b1e)), (counting(O40b2e) - counting_harsh(O40b2e)), (counting(O40b3e) - counting_harsh(O40b3e)), (counting(O50b1e) - counting_harsh(O50b1e)), (counting(O50b2e) - counting_harsh(O50b2e)), (counting(V1b1e) - counting_harsh(V1b1e)), (counting(V1b2e) - counting_harsh(V1b2e)), (counting(V40b1e) - counting_harsh(V40b1e)), (counting(V40b2e) - counting_harsh(V40b2e)), (counting(V40b3e) - counting_harsh(V40b3e)), (counting(V50b1e) - counting_harsh(V50b1e)), (counting(V50b2e) - counting_harsh(V50b2e)), (counting(V60b1e) - counting_harsh(V60b1e)), (counting(V60b2e) - counting_harsh(V60b2e)))

F1be = (counting(F1b1e)+counting(F1b2e))/2
F30be = (counting(F30b1e)+counting(F30b2e))/2
F40be = (counting(F40b1e)+counting(F40b2e)+counting(F40b3e))/3
O1be = (counting(O1b1e)+counting(O1b2e))/2
O30be = (counting(O30b1e)+counting(O30b2e))/2
O40be = (counting(O40b1e)+counting(O40b2e)+counting(O40b3e))/3
O50be = (counting(O50b1e)+counting(O50b2e))/2
V1be = (counting(V1b1e)+counting(V1b2e))/2
V40be = (counting(V40b1e)+counting(V40b2e)+counting(V40b3e))/3
V50be = (counting(V50b1e)+counting(V50b2e))/2
V60be = (counting(V60b1e)+counting(V60b2e))/2

F1bhe = (counting_harsh(F1b1e)+counting_harsh(F1b2e))/2
F30bhe = (counting_harsh(F30b1e)+counting_harsh(F30b2e))/2
F40bhe = (counting_harsh(F40b1e)+counting_harsh(F40b2e)+counting_harsh(F40b3e))/3
O1bhe = (counting_harsh(O1b1e)+counting_harsh(O1b2e))/2
O30bhe = (counting_harsh(O30b1e)+counting_harsh(O30b2e))/2
O40bhe = (counting_harsh(O40b1e)+counting_harsh(O40b2e)+counting_harsh(O40b3e))/3
O50bhe = (counting_harsh(O50b1e)+counting_harsh(O50b2e))/2
V1bhe = (counting_harsh(V1b1e)+counting_harsh(V1b2e))/2
V40bhe = (counting_harsh(V40b1e)+counting_harsh(V40b2e)+counting_harsh(V40b3e))/3
V50bhe = (counting_harsh(V50b1e)+counting_harsh(V50b2e))/2
V60bhe = (counting_harsh(V60b1e)+counting_harsh(V60b2e))/2

F1bde = F1be - F1bhe
F30bde = F30be - F30bhe
F40bde = F40be - F40bhe
O1bde = O1be - O1bhe
O30bde = O30be - O30bhe
O40bde = O40be - O40bhe
O50bde = O50be - O50bhe
V1bde = V1be - V1bhe
V40bde = V40be - V40bhe
V50bde = V50be - V50bhe
V60bde = V60be - V60bhe

epercents_harsh = (F1bhe, F30bhe, F40bhe, O1bhe, O30bhe, O40bhe, O50bhe, V1bhe, V40bhe, V50bhe, V60bhe)
epercents_dff = (F1bde, F30bde, F40bde, O1bde, O30bde, O40bde, O50bde, V1bde, V40bde, V50bde, V60bde)

####### BOTH

probs = gmm.predict_proba(FZ_tot)

probFZ = np.concatenate((FZ_tot, probs), axis = 1)

prob = Gmm.predict_proba(FE_tot)

probFE = np.concatenate((FE_tot, prob), axis = 1)

probabilities = np.concatenate((probs, prob), axis = 1)

F1b1 = probabilities[0:202,:]
F1b2 = probabilities[202:404,:]
F30b1 = probabilities[404:606,:]
F30b2 = probabilities[606:808,:]
F40b1 = probabilities[808:1010,:]
F40b2 = probabilities[1010:1212,:]
F40b3 = probabilities[1212:1414,:]
O1b1 = probabilities[1414:1616,:]
O1b2 = probabilities[1616:1818,:]
O30b1 = probabilities[1818:2020,:]
O30b2 = probabilities[2020:2222,:]
O40b1 = probabilities[2222:2424,:]
O40b2 = probabilities[2424:2626,:]
O40b3 = probabilities[2626:2828,:]
O50b1 = probabilities[2828:3030,:]
O50b2 = probabilities[3030:3232,:]
V1b1 = probabilities[3232:3434,:]
V1b2 = probabilities[3434:3636,:]
V40b1 = probabilities[3636:3838,:]
V40b2 = probabilities[3838:4040,:]
V40b3 = probabilities[4040:4242,:]
V50b1 = probabilities[4242:4444,:]
V50b2 = probabilities[4444:4646,:]
V60b1 = probabilities[4646:4848,:]
V60b2 = probabilities[4848:5050,:]

percentages_harsh = ((counting_harsh(F1b1)), (counting_harsh(F1b2)), (counting_harsh(F30b1)), (counting_harsh(F30b2)), (counting_harsh(F40b1)), (counting_harsh(F40b2)), (counting_harsh(F40b3)), (counting_harsh(O1b1)), (counting_harsh(O1b2)), (counting_harsh(O30b1)), (counting_harsh(O30b2)), (counting_harsh(O40b1)), (counting_harsh(O40b2)), (counting_harsh(O40b3)), (counting_harsh(O50b1)), (counting_harsh(O50b2)), (counting_harsh(V1b1)), (counting_harsh(V1b2)), (counting_harsh(V40b1)), (counting_harsh(V40b2)), (counting_harsh(V40b3)), (counting_harsh(V50b1)), (counting_harsh(V50b2)), (counting_harsh(V60b1)), (counting_harsh(V60b2)))

percentages_diff = ((counting(F1b1) - counting_harsh(F1b1)), (counting(F1b2) - counting_harsh(F1b2)), (counting(F30b1) - counting_harsh(F30b1)), (counting(F30b2) - counting_harsh(F30b2)), (counting(F40b1) - counting_harsh(F40b1)), (counting(F40b2) - counting_harsh(F40b2)), (counting(F40b3) - counting_harsh(F40b3)), (counting(O1b1) - counting_harsh(O1b1)), (counting(O1b2) - counting_harsh(O1b2)), (counting(O30b1) - counting_harsh(O30b1)), (counting(O30b2) - counting_harsh(O30b2)), (counting(O40b1) - counting_harsh(O40b1)), (counting(O40b2) - counting_harsh(O40b2)), (counting(O40b3) - counting_harsh(O40b3)), (counting(O50b1) - counting_harsh(O50b1)), (counting(O50b2) - counting_harsh(O50b2)), (counting(V1b1) - counting_harsh(V1b1)), (counting(V1b2) - counting_harsh(V1b2)), (counting(V40b1) - counting_harsh(V40b1)), (counting(V40b2) - counting_harsh(V40b2)), (counting(V40b3) - counting_harsh(V40b3)), (counting(V50b1) - counting_harsh(V50b1)), (counting(V50b2) - counting_harsh(V50b2)), (counting(V60b1) - counting_harsh(V60b1)), (counting(V60b2) - counting_harsh(V60b2)))

F1b = (counting(F1b1)+counting(F1b2))/2
F30b = (counting(F30b1)+counting(F30b2))/2
F40b = (counting(F40b1)+counting(F40b2)+counting(F40b3))/3
O1b = (counting(O1b1)+counting(O1b2))/2
O30b = (counting(O30b1)+counting(O30b2))/2
O40b = (counting(O40b1)+counting(O40b2)+counting(O40b3))/3
O50b = (counting(O50b1)+counting(O50b2))/2
V1b = (counting(V1b1)+counting(V1b2))/2
V40b = (counting(V40b1)+counting(V40b2)+counting(V40b3))/3
V50b = (counting(V50b1)+counting(V50b2))/2
V60b = (counting(V60b1)+counting(V60b2))/2
P1b = (counting(P1b1)+counting(P1b2))/2
P50b = (counting(P50b1)+counting(P50b2))/2

F1bh = (counting_harsh(F1b1)+counting_harsh(F1b2))/2
F30bh = (counting_harsh(F30b1)+counting_harsh(F30b2))/2
F40bh = (counting_harsh(F40b1)+counting_harsh(F40b2)+counting_harsh(F40b3))/3
O1bh = (counting_harsh(O1b1)+counting_harsh(O1b2))/2
O30bh = (counting_harsh(O30b1)+counting_harsh(O30b2))/2
O40bh = (counting_harsh(O40b1)+counting_harsh(O40b2)+counting_harsh(O40b3))/3
O50bh = (counting_harsh(O50b1)+counting_harsh(O50b2))/2
V1bh = (counting_harsh(V1b1)+counting_harsh(V1b2))/2
V40bh = (counting_harsh(V40b1)+counting_harsh(V40b2)+counting_harsh(V40b3))/3
V50bh = (counting_harsh(V50b1)+counting_harsh(V50b2))/2
V60bh = (counting_harsh(V60b1)+counting_harsh(V60b2))/2
P1bh = (counting_harsh(P1b1)+counting_harsh(P1b2))/2
P50bh = (counting_harsh(P50b1)+counting_harsh(P50b2))/2

F1bd = F1b - F1bh
F30bd = F30b - F30bh
F40bd = F40b - F40bh
O1bd = O1b - O1bh
O30bd = O30b - O30bh
O40bd = O40b - O40bh
O50bd = O50b - O50bh
V1bd = V1b - V1bh
V40bd = V40b - V40bh
V50bd = V50b - V50bh
V60bd = V60b - V60bh
P1bd = P1b - P1bh
P50bd = P50b - P50bh

percents_harsh = (F1bh, F30bh, F40bh, O1bh, O30bh, O40bh, O50bh, V1bh, V40bh, V50bh, V60bh)
percents_dff = (F1bd, F30bd, F40bd, O1bd, O30bd, O40bd, O50bd, V1bd, V40bd, V50bd, V60bd)

#######

index = np.arange(25)
index2 = np.arange(11)

labels = ('F1b1', 'F1b2', 'F30b1', 'F30b2', 'F40b1', 'F40b2', 'F40b3', 'O1b1', 'O1b2', 'O30b1', 'O30b2', 'O40b1', 'O40b2', 'O40b3', 'O50b1', 'O50b2', 'V1b1', 'V1b2', 'V40b1', 'V40b2', 'V40b3', 'V50b1', 'V50b2', 'V60b1', 'V60b2')
labelling = ('F1b', 'F30b', 'F40b', 'O1b', 'O30b', 'O40b', 'O50b', 'V1b', 'V40b', 'V50b', 'V60b')

##########
# PLOTTING
##########

fig = plt.figure(figsize = (8,8))
plt.scatter(FZ_tot[:,0], FZ_tot[:,1], c = y_k_means, s = 20, cmap = 'Wistia')
plt.scatter(ZAup, FAup, label = '4BW5 (Up)', color = 'white', marker = 'D', s = 60, edgecolor = 'black')
plt.scatter(ZAdn, FAdn, label = '4XDJ (Down)', color = 'grey', marker = 'D', s = 60, edgecolor = 'black')
center = k_means.cluster_centers_
plt.scatter(center[:,0], center[:,1], c = 'black', s = 200, alpha = 0.5)
v = [3,18,3,16]
plt.axis(v)
plt.axhline(5.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axvline(11.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.xlabel('Zipper ($\AA$)', fontsize = '18')
plt.ylabel('Fenestration ($\AA$)', fontsize = '18')
plt.savefig("clustering-fz-kmeans-nopop.png", format='png', dpi=300)
plt.savefig("clustering-fz-kmeans-nopop.svg", format='svg', dpi=300)

plt.clf()

plt.figure(figsize = (8,8))
plot_gmm(gmm, FZ_tot)
plt.scatter(FZ_tot[:,0], FZ_tot[:,1], c=label, s=20, cmap = 'rainbow', alpha = 0.3)
plt.scatter(ZAup, FAup, label = '4BW5 (Up)', color = 'white', marker = 'D', s = 60, edgecolor = 'black')
plt.scatter(ZAdn, FAdn, label = '4XDJ (Down)', color = 'grey', marker = 'D', s = 60, edgecolor = 'black')
plt.axhline(5.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axvline(11.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axhline(6.5, linestyle = '--', color = 'orange')
plt.axvline(10.5, linestyle = '--', color = 'orange')
plt.axhline(6.1, linestyle = '--', color = 'green')
plt.axvline(11.2, linestyle = '--', color = 'green')
v = [3,18,3,16]
plt.axis(v)
plt.xlabel('Zipper ($\AA$)', fontsize = '18')
plt.ylabel('Fenestration ($\AA$)', fontsize = '18')
plt.savefig("clustering-fz-gmm-nopop.png", format='png', dpi=300)
plt.savefig("clustering-fz-gmm-nopop.svg", format='svg', dpi=300)

plt.clf()

fig = plt.figure(figsize = (8,8))
plt.scatter(FE_tot[:,0], FE_tot[:,1], c = y_k_means, s = 20, cmap = 'Wistia')
plt.scatter(EAup, FAup, label = '4BW5 (Up)', color = 'white', marker = 'D', s = 60, edgecolor = 'black')
plt.scatter(EAdn, FAdn, label = '4XDJ (Down)', color = 'grey', marker = 'D', s = 60, edgecolor = 'black')
center = k_means.cluster_centers_
plt.scatter(center[:,0], center[:,1], c = 'black', s = 200, alpha = 0.5)
v = [3,18,3,16]
plt.axis(v)
plt.axhline(5.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axvline(11.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.xlabel('Expansion ($\AA$)', fontsize = '18')
plt.ylabel('Fenestration ($\AA$)', fontsize = '18')
plt.savefig("clustering-fe-kmeans-nopop.png", format='png', dpi=300)
plt.savefig("clustering-fe-kmeans-nopop.svg", format='svg', dpi=300)

plt.clf()

plt.figure(figsize = (8,8))
plot_gmm(Gmm, FE_tot)
plt.scatter(FE_tot[:,0], FE_tot[:,1], c=label, s=20, cmap = 'rainbow', alpha = 0.3)
plt.scatter(EAup, FAup, label = '4BW5 (Up)', color = 'white', marker = 'D', s = 60, edgecolor = 'black')
plt.scatter(EAdn, FAdn, label = '4XDJ (Down)', color = 'grey', marker = 'D', s = 60, edgecolor = 'black')
plt.axhline(5.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axvline(11.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axhline(6.5, linestyle = '--', color = 'orange')
plt.axvline(10.5, linestyle = '--', color = 'orange')
plt.axhline(6.1, linestyle = '--', color = 'green')
plt.axvline(11.2, linestyle = '--', color = 'green')
v = [3,18,3,16]
plt.axis(v)
plt.xlabel('Expansion ($\AA$)', fontsize = '18')
plt.ylabel('Fenestration ($\AA$)', fontsize = '18')
plt.savefig("clustering-fe-gmm-nopop.png", format='png', dpi=300)
plt.savefig("clustering-fe-gmm-nopop.svg", format='svg', dpi=300)

plt.clf()

plt.figure(figsize = (35,5))
width = 0.4
plt.figure(figsize = (30,5))
plt.bar(index - width/2, zpercentages_harsh, width, tick_label = labels, color = 'green', label = 'Zipper/Fenestration')
plt.bar(index - width/2, zpercentages_diff, width, bottom = zpercentages_harsh, color = 'green', alpha = 0.5)
plt.bar(index + width/2, epercentages_harsh, width, tick_label = labels, color = 'red', label = 'Expansion/Fenestration')
plt.bar(index + width/2, epercentages_diff, width, bottom = epercentages_harsh, color = 'red', alpha = 0.5)
plt.ylabel('Percentage of simulation time in Up state')
plt.legend(loc='best')
plt.savefig("clusterbarsall-nopop.png", format='png', dpi=300)
plt.savefig("clusterbarsall-nopop.svg", format='svg', dpi=300)

plt.clf()

plt.figure(figsize = (15,5))
plt.bar(index2 - width/2, zpercents_harsh, width, tick_label = labelling, color = 'green', label = 'Zipper/Fenestration')
plt.bar(index2 - width/2, zpercents_dff, width, bottom = zpercents_harsh, color = 'green', alpha = 0.5)
plt.bar(index2 + width/2, epercents_harsh, width, tick_label = labelling, color = 'red', label = 'Expansion/Fenestration')
plt.bar(index2 + width/2, epercents_dff, width, bottom = epercents_harsh, color = 'red', alpha = 0.5)
plt.ylabel('Percentage of simulation time in Up state')
plt.legend(loc='best', fontsize = 'medium')
plt.savefig("clusterbarsavs-nopop.png", format='png', dpi=300)
plt.savefig("clusterbarsavs-nopop.svg", format='svg', dpi=300)
