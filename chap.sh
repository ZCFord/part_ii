#!/bin/bash

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
python ~/Documents/PART_II/USEFUL/scripts/lipid-contacts-coli-at.py -f md-c.gro -t md-c.xtc -l $5 -o $1
perl ~/Documents/PART_II/USEFUL/scripts/make_b_factor.pl md-c.pdb $1_$5-contacts.csv 1 $1-bfacts.pdb
bash ~/Documents/PART_II/USEFUL/scripts/lca_writer.sh $1_$5-contacts.csv
