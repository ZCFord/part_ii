#!/bin/bash

gmx make_ndx -f md-c.gro -o rdf <<EOD

a N | a C13 | a C14 | a C15 | a C12 | a C11 | a P | a O11 | a O12 | a O13 | a O14
a N | a C13 | a C14 | a C15 | a C12 | a C11 | a P | a O11 | a O12 | a O13 | a O14 | a H13A | a H13B | a H13C | a H14A | a H14B | a H14C | a H15A | a H15B | a H15C | a H12A | a H12B | a H11A | a H11B
13 &! 25
name 24 Headgroups
name 26 Tails
q
EOD

gmx rdf -f md-c.xtc -s ../md.tpr -ref 2 -sel 24 -n rdf.ndx -o headrdf -rmax 1.5

gmx rdf -f md-c.xtc -s ../md.tpr -ref 2 -sel 26 -n rdf.ndx -o tailrdf -rmax 1.5

