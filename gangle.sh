#!/bin/bash

gmx make_ndx -f mdord.pdb -o TM-sel.ndx <<EOD
3 & r 194 206 & chain A
3 & r 194 206 & chain B
3 & r 194 222 & chain A
3 & r 194 222 & chain B
3 & r 232 248 & chain A
3 & r 232 248 & chain B
3 & r 232 264 & chain A
3 & r 232 264 & chain B
3 & r 301 312 & chain A
3 & r 301 312 & chain B
3 & r 301 328 & chain A
3 & r 301 328 & chain B
q
EOD

gmx gangle -f md-c.xtc -s ../md.tpr -n TM-sel.ndx -tu ns -g1 vector -g2 vector -oall 2a.xvg -group1 24 -group2 26
gmx gangle -f md-c.xtc -s ../md.tpr -n TM-sel.ndx -tu ns -g1 vector -g2 vector -oall 2b.xvg -group1 25 -group2 27

gmx gangle -f md-c.xtc -s ../md.tpr -n TM-sel.ndx -tu ns -g1 vector -g2 vector -oall 3a.xvg -group1 28 -group2 30
gmx gangle -f md-c.xtc -s ../md.tpr -n TM-sel.ndx -tu ns -g1 vector -g2 vector -oall 3b.xvg -group1 29 -group2 31

gmx gangle -f md-c.xtc -s ../md.tpr -n TM-sel.ndx -tu ns -g1 vector -g2 vector -oall 4a.xvg -group1 32 -group2 34
gmx gangle -f md-c.xtc -s ../md.tpr -n TM-sel.ndx -tu ns -g1 vector -g2 vector -oall 4b.xvg -group1 33 -group2 35

