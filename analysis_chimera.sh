#!/bin/bash

# script for running analysis for chimera!
# start from analysis folder
# USAGE: 1 = structure/pressure

python ~/Documents/PART_II/USEFUL/scripts/MDAnalysis_bilayer-thickness.py -f md-c.gro -t md-c.xtc -o $1

python ~/Documents/PART_II/USEFUL/scripts/MDAnalysis_angles.py -f mdord.pdb -t md-c.xtc -o $1

python ~/Documents/PART_II/USEFUL/scripts/MDAnalysis_distance.py -f mdord.pdb -t md-c.xtc -o $1 -u ~/Documents/PART_II/USEFUL/structures/UP_correctpro.pdb -d ~/Documents/PART_II/USEFUL/structures/DN_correctH.pdb

python ~/Documents/PART_II/USEFUL/scripts/MDAnalysis_rmsd.py -f mdord.pdb -ft md-c.xtc -s ~/Documents/PART_II/USEFUL/structures/UP_correctpro.pdb -o $1

bash ~/Documents/PART_II/USEFUL/scripts/fatslim_apl_at.sh

python ~/Documents/PART_II/USEFUL/scripts/plot_apl.py -o $1

bash ~/Documents/PART_II/USEFUL/scripts/rdf.sh

python ~/Documents/PART_II/USEFUL/scripts/plot_rdf.py -o $1

bash ~/Documents/PART_II/USEFUL/scripts/order_parameter_popc.sh

python ~/Documents/PART_II/USEFUL/scripts/plot_order.py -n 18 -d Yes -o $1

bash ~/Documents/PART_II/USEFUL/scripts/gangle.sh

python ~/Documents/PART_II/USEFUL/scripts/gangle-plot.py -o $1

mkdir chap
cd chap
mkdir ions
mkdir water
mkdir lipid

gmx trjconv -f ../md-c.xtc -s ../../md.tpr -o md-chap.xtc -fit rot+trans <<EOD
1
0
EOD

cd ions
chap -f ../md-chap.xtc -s ../../../md.tpr -sel-pathway 1 -sel-solvent 14 -de-bandwidth 0.14 -hydrophob-bandwidth 0.45 -pf-chan-dir-vec 0 0 -1

python ~/Documents/PART_II/USEFUL/scripts/chapcompound.py -o 40b2 -sol ions

cd ../water
chap -f ../md-chap.xtc -s ../../../md.tpr -sel-pathway 1 -sel-solvent 16 -de-bandwidth 0.14 -hydrophob-bandwidth 0.45 -pf-chan-dir-vec 0 0 -1

cd ../../
python ~/Documents/PART_II/USEFUL/scripts/lipid-contacts-coli-at.py -f md-c.gro -t md-c.xtc -l POPC -o $1
perl ~/Documents/PART_II/USEFUL/scripts/make_b_factor.pl md-c.pdb $1_POPC-contacts.csv 1 $1-bfacts.pdb
bash ~/Documents/PART_II/USEFUL/scripts/lca_writer.sh $1_POPC-contacts.csv
