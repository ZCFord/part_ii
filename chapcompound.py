#!/usr/bin/env python2

# this version currently being modified by Zoe! 
#
#
# CHAP - The Channel Annotation Package
# 
# Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
# Stephen J. Tucker
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


################################################################################
# CONFIGURATION
################################################################################

# load libraries:
import json                             # read in JSON files
import numpy as np                      # manipulate numeric vectors
from matplotlib import pyplot as pl     # plotting facilities
import argparse                         # parse command line arguments
from matplotlib import rc, rcParams
from matplotlib.gridspec import GridSpec
import pandas
import math

rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'

# get parameters from user input:
parser = argparse.ArgumentParser()
parser.add_argument(
    "-filename",
    nargs = "?",
    const = "output.json",
    default = "output.json")

parser.add_argument("-dpi",
    nargs = "?",
    const = 1200,
    default = 1200,
    type = int)

parser.add_argument("-o",
    default = 'chap'
)

parser.add_argument("-sol",
    default = 'water'
)

args = parser.parse_args()

################################################################################
# DATA READ-IN
################################################################################

# load output data from JSON file:
with open(args.filename) as data_file:
    data = json.load(data_file)

################################################################################
# PANDAS
################################################################################

df = pandas.DataFrame(data["residueSummary"]["s"])

df["z_mean"] = data["residueSummary"]["z"]["mean"]
df["pore_mean"] = data["residueSummary"]["poreRadius"]["mean"]
df["x_mean"] = data["residueSummary"]["x"]["mean"]
df["y_mean"] = data["residueSummary"]["y"]["mean"]
df["name"] = data["residueSummary"]["name"]
df["chain"] = data["residueSummary"]["chain"]
df["resID"] = data["residueSummary"]["id"]

del df["min"]
del df["max"]
del df["var"]
del df["sd"]

df.columns = [u's_mean', u'z_mean', u'pore_mean', u'x_mean', u'y_mean', u'name', u'chain', u'resID']

df["resID"] = df["resID"] + 73
chain_a = df[["s_mean", "z_mean", "pore_mean", "x_mean", "y_mean", "name", "resID", "chain"]][:261]
chain_b = df[["s_mean", "z_mean", "pore_mean", "x_mean", "y_mean", "name", "resID", "chain"]][261:]
chain_a["chain"] = "A"
chain_b["chain"] = "B"
chain_b["resID"] = chain_b["resID"] - 261
full = pandas.concat([chain_a, chain_b])

sf1 = full[["s_mean", "z_mean", "name", "resID", "chain"]][99:104]
sf2 = full[["s_mean", "z_mean", "name", "resID", "chain"]][208:213]
sf3 = full[["s_mean", "z_mean", "name", "resID", "chain"]][360:365]
sf4 = full[["s_mean", "z_mean", "name", "resID", "chain"]][469:474]
sf = pandas.concat([sf1, sf2, sf3, sf4])
thr = sf[sf["name"] == "THR"]
thr.groupby(["name"]).mean()
thr_av = thr.groupby(["name"]).mean()
gly = sf[sf["name"] == "GLY"]
gly1 = gly[gly["resID"] == 176]
gly2 = gly[gly["resID"] == 285]
gly = pandas.concat([gly1, gly2])
gly_av = gly.groupby(["name"]).mean()
avs = pandas.concat([gly_av, thr_av])
avs = avs[["s_mean", "z_mean"]]

full["x_2"] = full["x_mean"] * full["x_mean"]
full["y_2"] = full["y_mean"] * full["y_mean"]
full["xy_sum"] = full["x_2"] + full["y_2"]
full["xy_av"] = full["xy_sum"]**0.5

################################################################################
# THE COMBINED PLOT!
################################################################################

gs = GridSpec(1,10)
fig = pl.figure(figsize=(24,8))

##########

Q = fig.add_subplot(gs[:,0:3])

pl.plot((full["x_mean"])*10, (full["s_mean"])*10)
x = np.arange(20,110,10)
pl.fill_between(x, (avs.at['THR', 's_mean'])*10, (avs.at['GLY', 's_mean'])*10, facecolor = 'red', alpha = 0.2)
z = [20,90,-50,40]
pl.axis(z)
pl.xlabel("x coordinate ($\AA$)", fontsize = 18)
pl.ylabel("Distance along pore ($\AA$)", fontsize = 18)
pl.title("Protein", fontsize = 24)
pl.setp(Q.get_xticklabels(), fontsize = 14)
pl.setp(Q.get_yticklabels(), fontsize = 14)

##########

R = fig.add_subplot(gs[:,3:5])

pl.plot(
    	(np.array(data["pathwayProfile"]["radiusMean"]))*10,
	(np.array(data["pathwayProfile"]["s"]))*10,
	"k-")

v = [0.0,10.0,-30,15]
w = [0.0, 2.5, 5.0, 7.5, 10.0]
pl.axis(v)
pl.xticks(w)
pl.title("Radius Profile", fontsize = 24)
pl.setp(R.get_xticklabels(), fontsize = 14)
pl.setp(R.get_yticklabels(), fontsize = 14)
pl.xlabel("Radius ($\AA$)", fontsize = 18)

##########

P = fig.add_subplot(gs[:,5:])

s = np.array(data["pathwayProfileTimeSeries"]["s"])*10
t = np.array(data["pathwayProfileTimeSeries"]["t"])
t2 = t/1000
n = np.array(data["pathwayProfileTimeSeries"]["density"])

num_t = np.size(np.unique(t2))
num_s = np.size(np.unique(s))

S = s.reshape(num_t, num_s)
T = t2.reshape(num_t, num_s)
N = n.reshape(num_t, num_s)

pl.pcolormesh(
    T,
    S,
    N,
    cmap = "Blues")

if args.sol == 'water' :
    pl.clim(0,40)
    cbar = pl.colorbar()
    cbar.ax.set_ylabel("Number of Water Molecules ($\mathrm{nm}^{-3}$)", fontsize = 18)

else :
    pl.clim(0,0.12)
    cbar = pl.colorbar()
    cbar.ax.set_ylabel("Number of K+ Ions ($\mathrm{nm}^{-3}$)", fontsize = 18)

pl.title("Number Density Profile over Time", fontsize = 24)
v=[0,100,-30,15]
pl.axis(v)
pl.xlabel("Time (ns)", fontsize = 18)
pl.setp(P.get_xticklabels(), fontsize = 14)
pl.setp(P.get_yticklabels(), fontsize = 14)

##########

pl.subplots_adjust(wspace = 0.6)

##########

pl.savefig("{}_{}_chap_compound.png".format(args.o, args.sol), format='png', dpi=300)
pl.savefig("{}_{}_chap_compound.svg".format(args.o, args.sol), format='svg', dpi=300)
