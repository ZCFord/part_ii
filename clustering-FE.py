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
	res_i = x.select_atoms(sel1 + ' and name CA').coordinates()
	res_j = x.select_atoms(sel2 + ' and name CA').coordinates()
	dis = np.linalg.norm(res_i - res_j)
	return dis

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

##### DFPC
### 1b

F1b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/analysis/md-c.xtc', F1b1FE)
F1b1FE = np.squeeze(F1b1FE)

### -30b
from matplotlib import rc, rcParams
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

##### POPC
### 1b

P1b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc', P1b1FE)
P1b1FE = np.squeeze(P1b1FE)

P1b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc', P1b2FE)
P1b2FE = np.squeeze(P1b2FE)

### -40b

P50b1FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', P50b1FE)
P50b1FE = np.squeeze(P50b1FE)

P50b2FE = []
F_E('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', P50b2FE)
P50b2FE = np.squeeze(P50b2FE)

##########
# CLUSTERING
##########

FE_tot = np.concatenate((F1b1FE, F30b1FE, F30b2FE, F40b1FE, F40b2FE, F40b3FE, O1b1FE, O1b2FE, O30b1FE, O30b2FE, O40b1FE, O40b2FE, O40b3FE, O50b1FE, V1b1FE, V1b2FE, V40b1FE, V40b2FE, V40b3FE, V50b1FE, V50b2FE, V60b1FE, V60b2FE, P1b1FE, P1b2FE, P50b1FE, P50b2FE), axis = 0)

k_means = KMeans(n_clusters = 3)
k_means.fit(FE_tot)
y_k_means = k_means.predict(FE_tot)

gmm = GMM(n_components = 2, covariance_type='full')
labels = gmm.fit(FE_tot).predict(FE_tot)

##########
# PLOTTING
##########

fig = plt.figure(figsize = (8,8))
plt.scatter(FE_tot[:,1], FE_tot[:,0], c = y_k_means, s = 20, cmap = 'Wistia')
plt.scatter(FAup, EAup, label = '4BW5 (Up)', color = 'white', marker = 'D', s = 60, edgecolor = 'black')
plt.scatter(FAdn, EAdn, label = '4XDJ (Down)', color = 'grey', marker = 'D', s = 60, edgecolor = 'black')
center = k_means.cluster_centers_
plt.scatter(center[:,1], center[:,0], c = 'black', s = 200, alpha = 0.5)
v = [3,18,3,16]
plt.axis(v)
plt.axhline(5.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axvline(11.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.ylabel('Expansion ($\AA$)', fontsize = '18')
plt.xlabel('Fenestration ($\AA$)', fontsize = '18')
plt.savefig("clustering-fe-kmeans.png", format='png', dpi=300)
plt.savefig("clustering-fe-kmeans.svg", format='svg', dpi=300)

plt.clf()

plt.figure(figsize = (8,8))
plot_gmm(gmm, FE_tot)
plt.scatter(FE_tot[:,0], FE_tot[:,1], c=label, s=20, cmap = 'rainbow', alpha = 0.3)
plt.scatter(FAup, EAup, label = '4BW5 (Up)', color = 'white', marker = 'D', s = 60, edgecolor = 'black')
plt.scatter(FAdn, EAdn, label = '4XDJ (Down)', color = 'grey', marker = 'D', s = 60, edgecolor = 'black')
plt.axhline(5.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axvline(11.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axhline(6.5, linestyle = '--', color = 'orange')
plt.axvline(10.5, linestyle = '--', color = 'orange')
plt.axhline(6.1, linestyle = '--', color = 'green')
plt.axvline(11.2, linestyle = '--', color = 'green')
v = [3,18,3,16]
plt.axis(v)
plt.ylabel('Expansion ($\AA$)', fontsize = '18')
plt.xlabel('Fenestration ($\AA$)', fontsize = '18')
plt.savefig("clustering-fe-gmm.png", format='png', dpi=300)
plt.savefig("clustering-fe-gmm.svg", format='svg', dpi=300)
