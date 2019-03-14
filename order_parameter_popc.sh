#!/bin/bash

gmx make_ndx -f md-c.gro -o order_index_1.ndx <<EOD

ri 522 - 835
24 & a C21
24 & a C22
24 & a C23
24 & a C24
24 & a C25
24 & a C26
24 & a C27
24 & a C28
24 & a C29
24 & a C210
24 & a C211
24 & a C212
24 & a C213
24 & a C214
24 & a C215
24 & a C216
24 & a C217
24 & a C218
del 0 - 23
del 0
q
EOD

gmx make_ndx -f md-c.gro -o upper_index_1.ndx <<EOD

ri 522 - 678
24 & a C21
24 & a C22
24 & a C23
24 & a C24
24 & a C25
24 & a C26
24 & a C27
24 & a C28
24 & a C29
24 & a C210
24 & a C211
24 & a C212
24 & a C213
24 & a C214
24 & a C215
24 & a C216
24 & a C217
24 & a C218
del 0 - 23
del 0
q
EOD

gmx make_ndx -f md-c.gro -o lower_index_1.ndx <<EOD

ri 679 - 835
24 & a C21
24 & a C22
24 & a C23
24 & a C24
24 & a C25
24 & a C26
24 & a C27
24 & a C28
24 & a C29
24 & a C210
24 & a C211
24 & a C212
24 & a C213
24 & a C214
24 & a C215
24 & a C216
24 & a C217
24 & a C218
del 0 - 23
del 0
q
EOD

gmx make_ndx -f md-c.gro -o order_index_2.ndx <<EOD

ri 522 - 835
24 & a C31
24 & a C32
24 & a C33
24 & a C34
24 & a C35
24 & a C36
24 & a C37
24 & a C38
24 & a C39
24 & a C310
24 & a C311
24 & a C312
24 & a C313
24 & a C314
24 & a C315
24 & a C316
del 0 - 23
del 0
q
EOD

gmx make_ndx -f md-c.gro -o upper_index_2.ndx <<EOD

ri 522 - 678
24 & a C31
24 & a C32
24 & a C33
24 & a C34
24 & a C35
24 & a C36
24 & a C37
24 & a C38
24 & a C39
24 & a C310
24 & a C311
24 & a C312
24 & a C313
24 & a C314
24 & a C315
24 & a C316
del 0 - 23
del 0
q
EOD

gmx make_ndx -f md-c.gro -o lower_index_2.ndx <<EOD

ri 679 - 835
24 & a C31
24 & a C32
24 & a C33
24 & a C34
24 & a C35
24 & a C36
24 & a C37
24 & a C38
24 & a C39
24 & a C310
24 & a C311
24 & a C312
24 & a C313
24 & a C314
24 & a C315
24 & a C316
del 0 - 23
del 0
q
EOD

gmx order -f md-c.xtc -s ../md.tpr -n order_index_1.ndx -od order_graph_1.xvg

gmx order -f md-c.xtc -s ../md.tpr -n upper_index_1.ndx -od upper_graph_1.xvg

gmx order -f md-c.xtc -s ../md.tpr -n lower_index_1.ndx -od lower_graph_1.xvg

gmx order -f md-c.xtc -s ../md.tpr -n order_index_2.ndx -od order_graph_2.xvg

gmx order -f md-c.xtc -s ../md.tpr -n upper_index_2.ndx -od upper_graph_2.xvg

gmx order -f md-c.xtc -s ../md.tpr -n lower_index_2.ndx -od lower_graph_2.xvg

echo Have you checked the residue IDs of the lipids? Have you checked the number of carbons? Have you located md.tpr?
