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

F1b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/DFPC/Backwards/md/1bar/analysis/md-c.xtc', F1b1FZ)
F1b1FZ = np.squeeze(F1b1FZ)

### -30b
from matplotlib import rc, rcParams
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

##### POPC
### 1b

P1b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/1bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/1bar/repeat1/analysis/md-c.xtc', P1b1FZ)
P1b1FZ = np.squeeze(P1b1FZ)

P1b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/1bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/1bar/repeat2/analysis/md-c.xtc', P1b2FZ)
P1b2FZ = np.squeeze(P1b2FZ)

### -40b

P50b1FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc', P50b1FZ)
P50b1FZ = np.squeeze(P50b1FZ)

P50b2FZ = []
F_Z('/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/50bar/repeat2/analysis/mdord.pdb', '/sansom/s150/pemb4066/Documents/PART_II/TREK2/CG_T2/POPC/Backwards/md/50bar/repeat2/analysis/md-c.xtc', P50b2FZ)
P50b2FZ = np.squeeze(P50b2FZ)

##########
# CLUSTERING
##########

FZ_tot = np.concatenate((F1b1FZ, F30b1FZ, F30b2FZ, F40b1FZ, F40b2FZ, F40b3FZ, O1b1FZ, O1b2FZ, O30b1FZ, O30b2FZ, O40b1FZ, O40b2FZ, O40b3FZ, O50b1FZ, V1b1FZ, V1b2FZ, V40b1FZ, V40b2FZ, V40b3FZ, V50b1FZ, V50b2FZ, V60b1FZ, V60b2FZ, P1b1FZ, P1b2FZ, P50b1FZ, P50b2FZ), axis = 0)

k_means = KMeans(n_clusters = 3)
k_means.fit(FZ_tot)
y_k_means = k_means.predict(FZ_tot)

gmm = GMM(n_components = 2, covariance_type='full')
labels = gmm.fit(FZ_tot).predict(FZ_tot)

##########
# PLOTTING
##########

fig = plt.figure(figsize = (8,8))
plt.scatter(FZ_tot[:,1], FZ_tot[:,0], c = y_k_means, s = 20, cmap = 'Wistia')
plt.scatter(FAup, ZAup, label = '4BW5 (Up)', color = 'white', marker = 'D', s = 60, edgecolor = 'black')
plt.scatter(FAdn, ZAdn, label = '4XDJ (Down)', color = 'grey', marker = 'D', s = 60, edgecolor = 'black')
center = k_means.cluster_centers_
plt.scatter(center[:,1], center[:,0], c = 'black', s = 200, alpha = 0.5)
v = [3,18,3,16]
plt.axis(v)
plt.axhline(5.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axvline(11.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.ylabel('Zipper ($\AA$)', fontsize = '18')
plt.xlabel('Fenestration ($\AA$)', fontsize = '18')
plt.savefig("clustering-fz-kmeans.png", format='png', dpi=300)
plt.savefig("clustering-fz-kmeans.svg", format='svg', dpi=300)

plt.clf()

plt.figure(figsize = (8,8))
plot_gmm(gmm, FZ_tot)
plt.scatter(FZ_tot[:,0], FZ_tot[:,1], c=label, s=20, cmap = 'rainbow', alpha = 0.3)
plt.scatter(FAup, ZAup, label = '4BW5 (Up)', color = 'white', marker = 'D', s = 60, edgecolor = 'black')
plt.scatter(FAdn, ZAdn, label = '4XDJ (Down)', color = 'grey', marker = 'D', s = 60, edgecolor = 'black')
plt.axhline(5.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axvline(11.5, linestyle = '--', color = 'black', alpha = 0.4)
plt.axhline(6.5, linestyle = '--', color = 'orange')
plt.axvline(10.5, linestyle = '--', color = 'orange')
plt.axhline(6.1, linestyle = '--', color = 'green')
plt.axvline(11.2, linestyle = '--', color = 'green')
v = [3,18,3,16]
plt.axis(v)
plt.ylabel('Zipper ($\AA$)', fontsize = '18')
plt.xlabel('Fenestration ($\AA$)', fontsize = '18')
plt.savefig("clustering-fz-gmm.png", format='png', dpi=300)
plt.savefig("clustering-fz-gmm.svg", format='svg', dpi=300)
