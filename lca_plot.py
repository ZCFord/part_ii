import MDAnalysis
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import argparse
import csv

rcParams['axes.labelsize'] = 24
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = False

parser = argparse.ArgumentParser('Plotting LCA')

parser.add_argument(
	"-a",
	type = str,
	help = "1bar1"
)

parser.add_argument(
	"-b",
	type = str,	
	help = "1bar2"
)

parser.add_argument(
	"-c",
	type = str,
	help = "Stretch1"
)

parser.add_argument(
	"-d",
	type = str,	
	help = "Stretch2"
)

parser.add_argument(
	"-e",
	type = str,	
	help = "Stretch3"
)

parser.add_argument(
	"-o",
	type = str,
	default = "LCA",
	help = "output name"
)

args = parser.parse_args()

x = []
y = []

with open(args.a) as file_a:
	File_A = csv.reader(file_a, delimiter=",")
	headers = next(File_A)
	for row in File_A:
		x.append(float(row[0]))
		y.append(float(row[1]))

x = np.array(x)

q = []

with open(args.b) as file_b:
	File_B = csv.reader(file_b, delimiter=",")
	headers = next(File_B)
	for row in File_B:
		q.append(float(row[1]))

X = []
Y = []

with open(args.c) as file_c:
	File_C = csv.reader(file_c, delimiter=",")
	headers = next(File_C)
	for row in File_C:
		X.append(float(row[0]))
		Y.append(float(row[1]))

X = np.array(X)

Q = []

with open(args.d) as file_d:
	File_D = csv.reader(file_d, delimiter=",")
	headers = next(File_D)
	for row in File_D:
		Q.append(float(row[1]))

P = []

#with open(args.e) as file_e:
#	File_E = csv.reader(file_e, delimiter=",")
#	headers = next(File_E)
#	for row in File_E:
#		P.append(float(row[1]))

yq = np.stack((y,q), axis=1)
yqmean = np.array(yq.mean(axis=1))

YQ = np.stack((Y,Q), axis=1)
YQmean = np.array(YQ.mean(axis=1))

#YQP = np.stack((Y,Q,P), axis=1)
#YQPmean = np.array(YQP.mean(axis=1))

dif = np.subtract(yqmean, YQmean)
dif = np.array(dif)

#dif = np.subtract(yqmean, YQPmean)
#dif = np.array(dif)

x1 = x[0:261]+72
x2 = x[261:]-189

xre = np.concatenate((x1, x2))

data = np.stack((xre, dif), axis = 1)

plt.figure(figsize = (20,16))
plt.suptitle('Lipid Contact Analysis', fontsize = '36')
A = plt.subplot(211)
plt.title('Chain A', fontsize = 30)
plt.axvspan(73, 99, facecolor = '#A020F0', alpha = 0.4)
plt.axvspan(184, 230, facecolor = '#0000FF', alpha = 0.4)
plt.axvspan(233, 265, facecolor = '#FFA500', alpha = 0.4)
plt.axvspan(300, 333, facecolor = '#008000', alpha = 0.4)
plt.axvspan(198, 199, facecolor = 'blue', alpha = 0.5)
plt.axvspan(324, 325, facecolor = 'blue', alpha = 0.5)
plt.axvspan(237, 238, facecolor = 'green', alpha = 0.5)
plt.axvspan(326, 327, facecolor = 'green', alpha = 0.5)
plt.axvspan(212, 213, facecolor = 'red', alpha = 0.5)
plt.axvspan(322, 323, facecolor = 'red', alpha = 0.5)
plt.bar(data[0:261, 0], data[0:261, 1], color = 'black')
v = [73,333,-1.0, 1.0]
plt.axis(v)

B = plt.subplot(212, sharex = A)
plt.title('Chain B', fontsize = 30)
plt.axvspan(73, 99, facecolor = '#A020F0', alpha = 0.4)
plt.axvspan(184, 230, facecolor = '#0000FF', alpha = 0.4)
plt.axvspan(233, 265, facecolor = '#FFA500', alpha = 0.4)
plt.axvspan(300, 333, facecolor = '#008000', alpha = 0.4)
plt.axvspan(198, 199, facecolor = 'blue', alpha = 0.5)
plt.axvspan(324, 325, facecolor = 'blue', alpha = 0.5)
plt.axvspan(237, 238, facecolor = 'green', alpha = 0.5)
plt.axvspan(326, 327, facecolor = 'green', alpha = 0.5)
plt.axvspan(212, 213, facecolor = 'red', alpha = 0.5)
plt.axvspan(322, 323, facecolor = 'red', alpha = 0.5)
plt.bar(data[261:, 0], data[261:, 1], color = 'black')
v = [73,333,-1.0, 1.0]
plt.axis(v)
plt.xlabel("Residue Number", fontsize = '24')
plt.savefig("{}_lca.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_lca.svg".format(args.o), format='svg', dpi=300)
