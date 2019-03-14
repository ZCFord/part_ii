#!/usr/bin/env python2

# this version currently being modified by Zoe! 
# run from chap folder
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
with open('ions/output.json') as data_file1:
    data1 = json.load(data_file1)

with open('water/output.json') as data_file2:
    data2 = json.load(data_file2)

################################################################################
# PANDAS
################################################################################

##### IONS

df1 = pandas.DataFrame(data1["residueSummary"]["s"])

df1["z_mean"] = data1["residueSummary"]["z"]["mean"]
df1["pore_mean"] = data1["residueSummary"]["poreRadius"]["mean"]
df1["x_mean"] = data1["residueSummary"]["x"]["mean"]
df1["y_mean"] = data1["residueSummary"]["y"]["mean"]
df1["name"] = data1["residueSummary"]["name"]
df1["chain"] = data1["residueSummary"]["chain"]
df1["resID"] = data1["residueSummary"]["id"]

del df1["min"]
del df1["max"]
del df1["var"]
del df1["sd"]

df1.columns = [u's_mean', u'z_mean', u'pore_mean', u'x_mean', u'y_mean', u'name', u'chain', u'resID']

df1["resID"] = df1["resID"] + 73
chain_a1 = df1[["s_mean", "z_mean", "pore_mean", "x_mean", "y_mean", "name", "resID", "chain"]][:261]
chain_b1 = df1[["s_mean", "z_mean", "pore_mean", "x_mean", "y_mean", "name", "resID", "chain"]][261:]
chain_a1["chain"] = "A"
chain_b1["chain"] = "B"
chain_b1["resID"] = chain_b1["resID"] - 261
full1 = pandas.concat([chain_a1, chain_b1])

sf11 = full1[["s_mean", "z_mean", "name", "resID", "chain"]][99:104]
sf21 = full1[["s_mean", "z_mean", "name", "resID", "chain"]][208:213]
sf31 = full1[["s_mean", "z_mean", "name", "resID", "chain"]][360:365]
sf41 = full1[["s_mean", "z_mean", "name", "resID", "chain"]][469:474]
sf1 = pandas.concat([sf11, sf21, sf31, sf41])
thr1 = sf1[sf1["name"] == "THR"]
thr1.groupby(["name"]).mean()
thr_av1 = thr1.groupby(["name"]).mean()
gly1 = sf1[sf1["name"] == "GLY"]
gly11 = gly1[gly1["resID"] == 176]
gly21 = gly1[gly1["resID"] == 285]
gly1 = pandas.concat([gly11, gly21])
gly_av1 = gly1.groupby(["name"]).mean()
avs1 = pandas.concat([gly_av1, thr_av1])
avs1 = avs1[["s_mean", "z_mean"]]

full1["x_2"] = full1["x_mean"] * full1["x_mean"]
full1["y_2"] = full1["y_mean"] * full1["y_mean"]
full1["xy_sum"] = full1["x_2"] + full1["y_2"]
full1["xy_av"] = full1["xy_sum"]**0.5

##### WATER

df2 = pandas.DataFrame(data2["residueSummary"]["s"])

df2["z_mean"] = data2["residueSummary"]["z"]["mean"]
df2["pore_mean"] = data2["residueSummary"]["poreRadius"]["mean"]
df2["x_mean"] = data2["residueSummary"]["x"]["mean"]
df2["y_mean"] = data2["residueSummary"]["y"]["mean"]
df2["name"] = data2["residueSummary"]["name"]
df2["chain"] = data2["residueSummary"]["chain"]
df2["resID"] = data2["residueSummary"]["id"]

del df2["min"]
del df2["max"]
del df2["var"]
del df2["sd"]

df2.columns = [u's_mean', u'z_mean', u'pore_mean', u'x_mean', u'y_mean', u'name', u'chain', u'resID']

df2["resID"] = df2["resID"] + 73
chain_a2 = df2[["s_mean", "z_mean", "pore_mean", "x_mean", "y_mean", "name", "resID", "chain"]][:261]
chain_b2 = df2[["s_mean", "z_mean", "pore_mean", "x_mean", "y_mean", "name", "resID", "chain"]][261:]
chain_a2["chain"] = "A"
chain_b2["chain"] = "B"
chain_b2["resID"] = chain_b2["resID"] - 261
full2 = pandas.concat([chain_a2, chain_b2])

sf12 = full2[["s_mean", "z_mean", "name", "resID", "chain"]][99:104]
sf22 = full2[["s_mean", "z_mean", "name", "resID", "chain"]][208:213]
sf32 = full2[["s_mean", "z_mean", "name", "resID", "chain"]][360:365]
sf42 = full2[["s_mean", "z_mean", "name", "resID", "chain"]][469:474]
sf2 = pandas.concat([sf12, sf22, sf32, sf42])
thr2 = sf2[sf2["name"] == "THR"]
thr2.groupby(["name"]).mean()
thr_av2 = thr2.groupby(["name"]).mean()
gly2 = sf2[sf2["name"] == "GLY"]
gly12 = gly2[gly2["resID"] == 176]
gly22 = gly2[gly2["resID"] == 285]
gly2 = pandas.concat([gly12, gly22])
gly_av2 = gly2.groupby(["name"]).mean()
avs2 = pandas.concat([gly_av2, thr_av2])
avs2 = avs2[["s_mean", "z_mean"]]

full2["x_2"] = full2["x_mean"] * full2["x_mean"]
full2["y_2"] = full2["y_mean"] * full2["y_mean"]
full2["xy_sum"] = full2["x_2"] + full2["y_2"]
full2["xy_av"] = full2["xy_sum"]**0.5

################################################################################
# THE COMBINED PLOT!
################################################################################

gs = GridSpec(1,10)
fig = pl.figure(figsize=(24,8))

##########

Q = fig.add_subplot(gs[:,0:3])

pl.plot((full1["x_mean"])*10, (full1["s_mean"])*10)
x = np.arange(20,110,10)
pl.fill_between(x, (avs1.at['THR', 's_mean'])*10, (avs1.at['GLY', 's_mean'])*10, facecolor = 'red', alpha = 0.2)
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
    	(np.array(data1["pathwayProfile"]["radiusMean"]))*10,
	(np.array(data1["pathwayProfile"]["s"]))*10,
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

##

s2 = np.array(data2["pathwayProfileTimeSeries"]["s"])*10
t2 = np.array(data2["pathwayProfileTimeSeries"]["t"])
t22 = t2/1000
n2 = np.array(data2["pathwayProfileTimeSeries"]["density"])

num_t2 = np.size(np.unique(t22))
num_s2 = np.size(np.unique(s2))

S2 = s2.reshape(num_t2, num_s2)
T2 = t22.reshape(num_t2, num_s2)
N2 = n2.reshape(num_t2, num_s2)

pl.pcolormesh(
    T2,
    S2,
    N2,
    cmap = "Blues", alpha = 1.0)

pl.clim(0,40)
cbar = pl.colorbar()
cbar.ax.set_ylabel("Number of Water Molecules ($\mathrm{nm}^{-3}$)", fontsize = 18)

##

s1 = np.array(data1["pathwayProfileTimeSeries"]["s"])*10
t1 = np.array(data1["pathwayProfileTimeSeries"]["t"])
t21 = t1/1000
n1 = np.array(data1["pathwayProfileTimeSeries"]["density"])

num_t1 = np.size(np.unique(t21))
num_s1 = np.size(np.unique(s1))

S1 = s1.reshape(num_t1, num_s1)
T1 = t21.reshape(num_t1, num_s1)
N1 = n1.reshape(num_t1, num_s1)

pl.pcolormesh(
    T1,
    S1,
    N1,
    cmap = "Reds")

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

pl.savefig("ions_water.png", format='png', dpi=300)
pl.savefig("ions_water.svg", format='svg', dpi=300)
