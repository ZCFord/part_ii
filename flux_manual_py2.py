#!/usr/bin/env python3

# Edited by Zoe Ford to be Python 2 compatible
# also added argparse to manually enter z boundaries - edited pore boundary definition section
#
# Copyright (c) 2018 Hannah Smith and Gianni Klesse
#
# Calculates the number flux of a given subset of solvent through a pore by
# splitting the pore into three domains - extracellular, pore and intracellular
# and summing the total extracellular -> pore -> intracellular events and the
# total intracellular -> pore -> extracellular events and subtracting one from
# the other. Final result is net transport of molecules in the positive z
# direction. z margins move both. z boundaries inwards towards the pore centre
# to exclude molecules around both entrances. Timestep should be chosen so that
# it is significantly smaller than the time each molecule spends in the pore.
# Monatomic molecules and polyatomic molecules are inputted seperately. An
# unlimited number of different monatomic molecules can be analysed with the
# monatomic_selection flag. One species of polyatomic molecule can also be
# analysed with or without monatomic molecules with the polyatomic_selection
# flag. Its position is tracked by the position of one of its constituent atoms,
# specified with the polyatomic_position flag. This should be an atom type that
# appears only once in each molecule and that is the heaviest atom in the
# molecule, for example the O in water. Output is csv file with timestamp (ps),
# transport rate (no molecules/timestep), total quantity transported up to that
# timestamp (no molecules).
#
# sample use for all water molecules:
# python ~/scripts/flux_manual.py -s production.gro -t production.xtc -polyatomic_selection 'resname SOL' -polyatomic_position 'type O' -o SOL -z_margin 10 -dt 10
# single ion:
# python ~/scripts/flux_manual.py -s production.gro -t production.xtc -monatomic_selection 'resid 14030' -o ID14030 -z_margin 10 -dt 10
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
parser = argparse.ArgumentParser('Calculates flux of given molecule across bilayer.')
parser.add_argument(
    '-s',
    type = str,
    default = 'production.gro',
    help = 'Structure file.'
)
parser.add_argument(
    '-t',
    type = str,
    default = 'production.trr',
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
    default = 'event_data',
    help = 'Output file.'
)
parser.add_argument(
    '-z_margin',
    type = float,
    default = 0,
    help = 'Safety margin on z limits in Angstroms.'
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

def solvent_domain_array_calc(sel, z_lo, z_hi, timestep, idx):
    """
    Calculate an array of solvent domain for a given timestamp.

    Keyword argument:
    sel -- set of molecules being analysed
    z_lo -- lower z boundary of protein
    z_hi -- higher z boundary of protein
    timestep -- time between frames (ps)
    """
    timestamp = u.trajectory.time
    solvent_z_coord = sel.positions[:, 2]
    domain_snapshot = np.array([assign_domain_index(z_mol, z_lo, z_hi) for z_mol in solvent_z_coord])
    solvent_domain_array[idx][0] = timestamp 
    solvent_domain_array[idx][1:] = domain_snapshot

def assign_domain_index(z_mol, z_lo, z_hi):
    """
    Assign a positional domain 1=> intracellular 2=> pore 3=> extracellular
 
    Keyword argument:
    z_mol -- z com position of current molecule being analysed
    z_lo -- lower z boundary of protein
    z_hi -- higher z boundary of protein
    """
    if z_mol < z_lo:
        return -1
    elif z_mol > z_hi:
        return 1
    else:
        return 0

################################################################################ 
# TRAJECTORY ANALYSIS
################################################################################ 
          
# creates a universe from given structure and trajectory files
u = mda.Universe(args.s, args.t)

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

# get protein boundaries in z-dir
if args.u == None and args.l == None :
    protein = u.select_atoms("protein")
    bbox = protein.bbox()
    z_lo = bbox[0][2] + args.z_margin
    z_hi = bbox[1][2] - args.z_margin
else :
    z_lo = args.l
    z_hi = args.u
print (z_lo)
print (z_hi)

# calculate frame parameters
# set default start and end to that of entire trajectory
if args.b == None:
    b = u.trajectory[0].time
else:
    b = args.b
if args.e == None:
    e = u.trajectory[-1].time
else:
    e = args.e
dt = args.dt
no_frames = int((e-b)/dt) + 1
start_frame = int(b/u.trajectory.dt)
end_frame = int(e/u.trajectory.dt)

# get number of mol in subset
no_solvent_mol = sel.n_residues

# initialises solvent domain array
solvent_domain_array = np.zeros((no_frames, no_solvent_mol+1), dtype=int)

# select each solvent molecule as individual atom group:
timestep = u.trajectory.dt

# runs func that writes domain data to output file
idx = 0
for ts in u.trajectory[start_frame:end_frame+1:int(dt/u.trajectory.dt)]:
    solvent_domain_array_calc(sel, z_lo, z_hi, timestep, idx)
    idx += 1

# write to CSV output:
np.savetxt(
    "{}_zm={}_b={}_domain_data.csv".format(args.o, args.z_margin, b),
    solvent_domain_array,
    comments = "",
    delimiter=",")

event_time_array = np.asarray([[0, 0]])

for mol in range(1, no_solvent_mol+1):
    # initialise events totals and mol's trajectory 
    total_up_events = 0
    total_down_events = 0
    net_events = 0
    single_mol_domain_array = np.empty([no_frames, 2])
    single_mol_domain_array[:,0] = solvent_domain_array[:,0]
    single_mol_domain_array[:,1] = solvent_domain_array[:,mol]
    
    # exclude start in pore
    while single_mol_domain_array[0,1] == 0:
        single_mol_domain_array = single_mol_domain_array[1:]
        if len(single_mol_domain_array == 0):
            break
    
    # calc positional changes
    domain_change_array = np.empty([len(single_mol_domain_array)-1, 2])
    domain_change_array[:,0] = single_mol_domain_array[1:,0]
    domain_change_array[:,1] = np.diff(single_mol_domain_array[:,1], axis=0)

    # remove zeros
    domain_change_nozeros = domain_change_array[np.nonzero(domain_change_array[:,1])]   
    # split trajectory at every periodic boundary jump
    split_domain_change_array = np.split(domain_change_nozeros, np.where(abs(domain_change_nozeros[:,1]) == 2)[0])
    for split_trajectory in split_domain_change_array:
        # ignore empty trajectories
        if len(split_trajectory[:,1]) < 2:
            continue

        # remove periodic boundary jumps
        split_trajectory = split_trajectory[abs(split_trajectory[:,1]) != 2]

        # make trajectory length even
        if len(split_trajectory[:,1])%2 == 1:
            split_trajectory = split_trajectory[:-1]
        
        # count total events in both directions
        events_array = np.empty([int(len(split_trajectory)/2), 2])
        events_array[:,0] = split_trajectory[1::2,0]
        events_array[:,1] = (split_trajectory[0::2,1] + split_trajectory[1::2,1])/2
        events_array_nozeros = events_array[np.nonzero(events_array[:,1])]
        if(events_array_nozeros.size != 0):
            event_time_array = np.append(event_time_array, events_array_nozeros, axis=0)

event_time_array = event_time_array[1:]

results = np.zeros([no_frames,3])
results[:,0] = solvent_domain_array[:,0]

for idx in range(0,len(event_time_array[:,0])):
    results[np.where(results[:,0]==event_time_array[idx, 0]), 1] += event_time_array[idx, 1]

results[:, 2] = np.asarray(np.cumsum(results[:, 1]))

# write to CSV output:
np.savetxt(
    "{}_zm={}_b={}_manual.csv".format(args.o, args.z_margin, b),
    results,
    comments = "",
    #header = "time, transport_rate, transport_quantity",
    delimiter=",")

