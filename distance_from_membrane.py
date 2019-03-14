import MDAnalysis
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis.analysis.leaflet

u = MDAnalysis.Universe('./md-c.pdb', './md-100.xtc')
L = MDAnalysis.analysis.leaflet.LeafletFinder(u, 'name P*')
downP = L.groups(0)

distanceA = []
distanceB = []
pcom = []

for ts in u.trajectory:
	alaA = u.select_atoms('resid 337 and segid A and name CA').coordinates()
	alaB = u.select_atoms('resid 337 and segid B and name CA').coordinates()
	downCOM = downP.center_of_mass()
	distanceA.append(downCOM[2] - alaA[0,2])
	distanceB.append(downCOM[2] - alaB[0,2])
	pcom.append(downCOM[2])

plt.plot(distanceA)
plt.plot(distanceB)
