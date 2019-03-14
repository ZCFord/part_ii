#!/bin/bash

##FATSLiM APL COARSE-GRAIN
## will only take .gro files

gmx make_ndx -f md500.gro -o apl.ndx <<EOD
a PO4
name 16 headgroups
q
EOD

mkdir apl_frames

fatslim apl -c md500.gro -n apl.ndx -t md500.xtc --plot-apl apl.xvg --export-apl-raw apl_frames/apl.csv -e 300000

xmgrace apl.xvg
