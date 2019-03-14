#!/bin/bash

python ~/Documents/PART_II/USEFUL/scripts/lipid-contacts-coli-at.py -f md-c.gro -t md-c.xtc -l $2 -o $1
perl ~/Documents/PART_II/USEFUL/scripts/make_b_factor.pl md-c.pdb $1_$2-contacts.csv 1 $1-bfacts.pdb
bash ~/Documents/PART_II/USEFUL/scripts/lca_writer.sh $1_$2-contacts.csv
