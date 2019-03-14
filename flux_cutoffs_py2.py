#!/usr/bin/env python3

# Edited by Zoe Ford to be Python 2 compatible
# also added argparse to manually enter z boundaries - edited pore boundary definition section
#
# Copyright (c) 2018 Hannah Smith and Gianni Klesse
# 
# Calculates the number flux of a given subset of solvent through a cylindrical
# pore directly from their velocities. Displacement of a molecule between
# subsequent frames is used as a pseudo velocity to smooth thermal noise.
# Transport rate is calculated from the sum of the displacements of molecules
# in the pore and transport quantity as the cumulative sum. z margins move both
# z boundaries inwards towards the pore centre to exclude molecules around both entrances.
# Pore radius is user input to make the cylindrical cutoff more biologically
# realistic. Timestep should be chosen so that it is significantly smaller than
# the time each molecule spends in the pore. Monatomic molecules and polyatomic
# molecules are inputted seperately. An unlimited number of different monatomic
# molecules can be analysed with the monatomic_selection flag. One species of
# polyatomic molecule can also be analysed with or without monatomic molecules
# with the polyatomic_selection flag. Its position is tracked by the position of
# one of its constituent atoms, specified with the polyatomic_position flag.
# This should be an atom type that appears only once in each molecule and that
# is the heaviest atom in the molecule, for example the O in water. Cuts off z
# displacement at pore boundaries so that when molecule is crossing z boundaries
# only displacement inside the pore contributes. Output is csv file with
# timestamp (ps), transport rate (no molecules/timestep), total quantity
# transported up to that timestamp (no molecules).
#
# sample use for all water molecules:
# python ~/scripts/flux_cutoffs.py -s production.gro -t production.xtc -polyatomic_selection 'resname SOL' -polyatomic_position 'type O' -o SOL -z_margin 10 -radius 10 -dt 10
# single ion:
# python ~/scripts/flux_cutoffs.py -s production.gro -t production.xtc -monatomic_selection 'resid 14030' -o ID14030 -z_margin 10 -radius 10 -dt 10
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


import MDAnalysis as mda
import numpy as np
from numpy.linalg import norm
import json
import argparse
import csv


###############################################################################
# ARGUMENT HANDLING
###############################################################################

# parse command line arguments:
parser = argparse.ArgumentParser('Calculates flux of given subset across'
    ' bilayer directly from velocities.')
parser.add_argument(
    '-s',
    type = str,
    default = 'production.gro',
    help = 'Structure file.'
)
parser.add_argument(
    '-t',
    type = str,
    default = 'production.xtc',
    help = 'Trajectory file.'
)
parser.add_argument(
    '-b',
    type = float,
    default = None,
    help = 'Start time in picoseconds.')
parser.add_argument(
    '-e',
    type = float,
    default = None,
    help = 'End time in picoseconds.')
parser.add_argument(
    '-dt',
    type = float,
    default = 100,
    help = 'Time step in picoseconds.')
parser.add_argument(
    '-monatomic_selection',
    type = str,
    default = None,
    help = 'Monatomic solvent being analysed.'
)
parser.add_argument(
    '-polyatomic_selection',
    type = str,
    default = None,
    help = 'Polyatomic solvent being analysed.'
)
parser.add_argument(
    '-polyatomic_position',
    type = str,
    default = 'type O',
    help = 'Type of atom in polyatomic selection whose position is tracked.'
)
parser.add_argument(
    '-o',
    type = str,
    default = 'velocities_data',
    help = 'Output file.'
)
parser.add_argument(
    '-z_margin',
    type = float,
    default = 0,
    help = 'Safety margin on z limits in Angstroms.'
)
parser.add_argument(
    '-radius',
    type = float,
    default = 5,
    help = 'Pore radius in Angstroms.'
)
parser.add_argument(
    '-u',
    type = float,
    default = None,
    help = 'Upper z boundary in Angstroms.'
)
parser.add_argument(
    '-l',
    type = float,
    default = None,
    help = 'Lower z boundary in Angstroms.'
)
args = parser.parse_args()


################################################################################ 
# FUNCTION DEFINTIONS
################################################################################ 

def is_outside_pore(coords, cog, pore_radius):
    """
    Calculate array of booleans specifying whether the radius of a molecule is
    in the pore.
    Dimensions of output: no frames x no residues.

    Keyword arguments:
    coords -- array of x, y, z coords for each mol and timestep in Angstroms
    cog -- center of geometry of pore in Angstroms
    pore_radius -- radius of pore in Angstroms
    """
    return ~(((coords[:,:,0]-cog[0])**2 + (coords[:,:,1]-cog[1])**2) < pore_radius**2)


################################################################################ 
# TRAJECTORY ANALYSIS
################################################################################ 

# create a universe from given structure and trajectory files
u = mda.Universe(args.s, args.t)

# pore boundary definition
if args.u == None and args.l == None :
    protein = u.select_atoms("protein")
    bbox = protein.bbox()
    cog = protein.center_of_geometry()
    z_lo = bbox[0][2] + args.z_margin
    z_hi = bbox[1][2] - args.z_margin
else :
    protein = u.select_atoms("protein")
    bbox = protein.bbox()
    cog = protein.center_of_geometry()
    z_lo = args.l
    z_hi = args.u    

pore_radius = args.radius

# select specified solvent subset and sets default to O positions in water
if args.monatomic_selection == None and args.polyatomic_selection == None:
    sel = u.select_atoms("resname SOL and type O")
elif args.monatomic_selection != None and args.polyatomic_selection == None:
    sel = u.select_atoms(args.monatomic_selection)
elif args.monatomic_selection == None and args.polyatomic_selection != None:
    sel = u.select_atoms("{} and {}".format(args.polyatomic_selection, args.polyatomic_position))
else:
    sel = u.select_atoms("{} or {} and {}".format(args.monatomic_selection, args.polyatomic_selection, args.polyatomic_position))

# check that resulting list of atoms only contains one from each molecule
if(len(np.unique(sel.resids)) != len(sel.resids)):
    raise Exception('Selection contains duplicate molecules. '
    'Check monatomic selection does not include polyatomic molecules.')

# calculate frame parameters
# set default start and end to that of entire trajectory
if args.b == None or args.b < u.trajectory[0].time:
    b = u.trajectory[0].time
else:
    b = args.b
if args.e == None  or args.e > u.trajectory[-1].time:
    e = u.trajectory[-1].time
else:
    e = args.e
dt = args.dt
no_frames = int((e-b)/dt) + 1
start_frame = int(b/u.trajectory.dt)
end_frame = int(e/u.trajectory.dt)

# loop over entire trajectory to give spatial and temporal coordinates
# for each mol and frame
idx = 0
coords = np.empty([no_frames, sel.n_residues, 3])
time_array = np.empty(no_frames) 
for ts in u.trajectory[start_frame:end_frame+1:int(dt/u.trajectory.dt)]:
    time_array[idx] = ts.time
    coords[idx, :] = sel.positions
    idx += 1

# set copy of z coords to boundary if z coords is outside of pore boundaries
z_coords_cutoff = np.empty([no_frames, sel.n_residues])
z_coords_cutoff = coords[:, :, 2]
z_coords_cutoff[z_coords_cutoff < z_lo] = z_lo
z_coords_cutoff[z_coords_cutoff > z_hi] = z_hi
print (z_hi)
print (z_lo)

# calculate z diplacement, removes periodic boundary jumps and sets to zero if
# radius is larger than the pore radius
z_displacements = np.empty([no_frames, sel.n_residues])
z_displacements[0, :] = 0
z_displacements[1:, :] = np.diff(z_coords_cutoff, axis=0)
z_displacements[abs(z_displacements) >= z_hi - z_lo] = 0
z_displacements[is_outside_pore(coords, cog, pore_radius)] = 0

# calculate total and cumulative sum of z displacement of all molecules
results = np.empty([no_frames, 3])
results[:, 0] = time_array
results[:, 1] = np.asarray(np.sum(z_displacements, axis=1))/(z_hi - z_lo)
results[:, 2] = np.asarray(np.cumsum(results[:, 1]))

# write to CSV output:
np.savetxt(
    "{}_zm={}_r={}.csv".format(args.o, args.z_margin, b),
    results,
    comments = "",
    #header = "time,transport_rate,transport_quantity",
    delimiter=",")

