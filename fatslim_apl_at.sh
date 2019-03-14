#!/bin/bash

##FATSLiM APL ATOMISTIC
## will only take .gro files

gmx make_ndx -f md-c.gro -o apl.ndx <<EOD
a NC3 C12 C13 C14 C15 C11 P O13 O14 O12 O11
name 24 headgroups
q
EOD

mkdir apl_frames

fatslim apl -c md-c.gro -n apl.ndx -t md-c.xtc --plot-apl apl.xvg --export-apl-raw apl_frames/apl.csv

