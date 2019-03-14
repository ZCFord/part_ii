import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import csv
import argparse

rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'

#### ARGUMENT HANDLING

parser = argparse.ArgumentParser('Plots flux')

parser.add_argument(
	"-f",
	type = str,
	required=True,
	help="CSV file input"
)

parser.add_argument(
	"-o",
	type = str,
	default = None,
	help = "output name"
)

args = parser.parse_args()

#### ACTUAL THING

x=[]
y=[]

with open(args.f) as csvfile:
    plots = csv.reader(csvfile, delimiter=",")
    for row in plots:
        x.append(float(row[0])/1000)
        y.append(float(row[-1]))

plt.figure(figsize=(10,8))
plt.plot(x, y, label = 'Ion Movement', color = 'pink', linewidth = 2)
plt.title("Movement of K ions through the pore", fontsize = 24)
plt.ylabel("Flux", fontsize = '16')
plt.xlabel("Time (ns)", fontsize = '16')
plt.savefig("{}_flux.png".format(args.o), format='png', dpi=300)
plt.savefig("{}_flux.svg".format(args.o), format='svg', dpi=300)

