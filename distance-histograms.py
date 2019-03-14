import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'

##########
# DEFINITIONS
##########

def dis (sel1, sel2, x) :
	res_i = x.select_atoms(sel1 + ' and name CA').positions
	res_j = x.select_atoms(sel2 + ' and name CA').positions
	dis = np.linalg.norm(res_i - res_j)
	return dis

##########
# CALCULATIONS
##########

u = MDAnalysis.Universe('../../../../USEFUL/structures/DN_correctResid.pdb')
protein = u.select_atoms('protein')

FF = dis('resid 226 and segid A', 'resid 226 and segid B', protein)
GG = dis('resid 324 and segid A', 'resid 324 and segid B', protein)
PG = dis('resid 198 and segid A', 'resid 324 and segid B', protein)
GP = dis('resid 324 and segid A', 'resid 198 and segid B', protein)

##### Unstretched

u1 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/analysis/md-c.xtc')
protein1 = u1.select_atoms('protein')

u2 = MDAnalysis.Universe('../../DFPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DFPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
protein2 = u2.select_atoms('protein')

u3 = MDAnalysis.Universe('../../DOPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/analysis/md-c.xtc')
protein3 = u3.select_atoms('protein')

u4 = MDAnalysis.Universe('../../DOPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DOPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
protein4 = u4.select_atoms('protein')

u5 = MDAnalysis.Universe('../../DVPC/Backwards/md/1bar/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/analysis/md-c.xtc')
protein5 = u5.select_atoms('protein')

u6 = MDAnalysis.Universe('../../DVPC/Backwards/md/1bar/repeat/analysis/mdord.pdb', '../../DVPC/Backwards/md/1bar/repeat/analysis/md-c.xtc')
protein6 = u6.select_atoms('protein')

F2F = []
G2G = []
P2G = []

for ts in u1.trajectory:
    F2F.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein1))
    G2G.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein1))
    P2G.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein1))
    P2G.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein1))

for ts in u2.trajectory:
    F2F.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein2))
    G2G.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein2))
    P2G.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein2))
    P2G.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein2))

for ts in u3.trajectory:
    F2F.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein3))
    G2G.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein3))
    P2G.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein3))
    P2G.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein3))

for ts in u4.trajectory:
    F2F.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein4))
    G2G.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein4))
    P2G.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein4))
    P2G.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein4))

for ts in u5.trajectory:
    F2F.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein5))
    G2G.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein5))
    P2G.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein5))
    P2G.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein5))

for ts in u6.trajectory:
    F2F.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein6))
    G2G.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein6))
    P2G.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein6))
    P2G.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein6))

##### Stretched

u7 = MDAnalysis.Universe('../../DFPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DFPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
protein7 = u7.select_atoms('protein')

u8 = MDAnalysis.Universe('../../DOPC/Backwards/md/40bar/repeat2/analysis/mdord.pdb', '../../DOPC/Backwards/md/40bar/repeat2/analysis/md-c.xtc')
protein8 = u8.select_atoms('protein')

u9 = MDAnalysis.Universe('../../DOPC/Backwards/md/50bar/repeat1/analysis/mdord.pdb', '../../DOPC/Backwards/md/50bar/repeat1/analysis/md-c.xtc')
protein9 = u9.select_atoms('protein')

u10 = MDAnalysis.Universe('../../DFPC/Backwards/md/30bar/repeat1/analysis/mdord.pdb', '../../DFPC/Backwards/md/30bar/repeat1/analysis/md-c.xtc')
protein10 = u10.select_atoms('protein')

u11 = MDAnalysis.Universe('../../DVPC/Backwards/md/40bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/40bar/repeat1/analysis/md-c.xtc')
protein11 = u11.select_atoms('protein')

u12 = MDAnalysis.Universe('../../DVPC/Backwards/md/60bar/repeat1/analysis/mdord.pdb', '../../DVPC/Backwards/md/60bar/repeat1/analysis/md-c.xtc')
protein12 = u12.select_atoms('protein')

f2f = []
g2g = []
p2g = []

for ts in u7.trajectory:
    f2f.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein7))
    g2g.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein7))
    p2g.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein7))
    p2g.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein7))

for ts in u8.trajectory:
    f2f.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein8))
    g2g.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein8))
    p2g.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein8))
    p2g.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein8))

for ts in u9.trajectory:
    f2f.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein9))
    g2g.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein9))
    p2g.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein9))
    p2g.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein9))

for ts in u10.trajectory:
    f2f.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein10))
    g2g.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein10))
    p2g.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein10))
    p2g.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein10))

for ts in u11.trajectory:
    f2f.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein11))
    g2g.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein11))
    p2g.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein11))
    p2g.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein11))

for ts in u12.trajectory:
    f2f.append(dis('resid 226 and segid A', 'resid 226 and segid B', protein12))
    g2g.append(dis('resid 324 and segid A', 'resid 324 and segid B', protein12))
    p2g.append(dis('resid 198 and segid A', 'resid 324 and segid B', protein12))
    p2g.append(dis('resid 324 and segid A', 'resid 198 and segid B', protein12))

##########
# PLOTTING
##########

plt.figure(figsize = (18,6))

F = plt.subplot(131)
plt.hist(F2F, bins=19, normed = True, color = 'black')
plt.hist(f2f, bins=19, normed = True, color = 'red')
plt.axvline(FF, color = 'black', linewidth = 4, alpha = 0.3)
plt.setp(F.get_xticklabels(), fontsize = 14)
plt.setp(F.get_yticklabels(), fontsize = 14)
plt.xlabel('F226 - F226 ($\AA$)', fontsize = '14')
plt.ylabel('Probability', fontsize = '14')
v = [50, 70, 0, 0.3]
plt.axis(v)

G = plt.subplot(132)
plt.hist(G2G, bins=19, normed = True, color = 'black')
plt.hist(g2g, bins=19, normed = True, color = 'red')
plt.axvline(GG, color = 'black', linewidth = 4, alpha = 0.3)
plt.setp(G.get_xticklabels(), fontsize = 14)
plt.setp(G.get_yticklabels(), fontsize = 14)
plt.xlabel('G324 - G324 ($\AA$)', fontsize = '14')
plt.ylabel('Probability', fontsize = '14')
v = [10, 30, 0, 0.4]
plt.axis(v)

P = plt.subplot(133)
plt.hist(P2G, bins=19, normed = True, color = 'black')
plt.hist(p2g, bins=19, normed = True, color = 'red')
plt.axvline(PG, color = 'black', linewidth = 4, alpha = 0.3)
plt.axvline(GP, color = 'black', linewidth = 4, alpha = 0.3)
plt.setp(P.get_xticklabels(), fontsize = 14)
plt.setp(P.get_yticklabels(), fontsize = 14)
plt.xlabel('P198 - G324 ($\AA$)', fontsize = '14')
plt.ylabel('Probability', fontsize = '14')
v = [0, 15, 0, 0.6]
plt.axis(v)

plt.savefig("disthist.png", format='png', dpi=300)
plt.savefig("disthist.svg", format='svg', dpi=300)
