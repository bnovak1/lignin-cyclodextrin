"""
Get atom numbers from plumed file for lignin head and tail, and cyclodextrin oxygen atoms in the secondary face. Then write a Tcl file to draw arrows in VMD from lignin head to tail, and from the cyclodextrin center to the secondary face oxygen atoms.
"""

import argparse

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--plumed", type=str, help="Plumed input file")
parser.add_argument("--tcl", type=str, help="Output Tcl file")
parser.add_argument("--arrow_length", type=float, help="Length of the arrow")
args = parser.parse_args()

# Read from the plumed file
with open(args.plumed, "r", encoding="utf-8") as f:
    lines = f.readlines()

for line in lines:
    if "ligninhead:" in line:
        atoms_head = np.array(line.split("=")[-1].strip().split(",")).astype(int) - 1
    if "lignintail:" in line:
        atoms_tail = np.array(line.split("=")[-1].strip().split(",")).astype(int) - 1
    if "bcdsecondary:" in line:
        atoms_bcd_2 = np.array(line.split("=")[-1].strip().split(",")).astype(int) - 1

# Write the Tcl file
TCL_STR = "# Arrow & head length\n"
TCL_STR += f"set arrow_length {args.arrow_length}\n\n"

TCL_STR += "# Atom selections\n"
TCL_STR += f'set ligninhead [atomselect top "index {" ".join(atoms_head.astype(str))}"]\n'
TCL_STR += f'set lignintail [atomselect top "index {" ".join(atoms_tail.astype(str))}"]\n'
TCL_STR += f'set bcdsecondary [atomselect top "index {" ".join(atoms_bcd_2.astype(str))}"]\n'
TCL_STR += 'set bcd [atomselect top "resname BCDC"]\n\n'

TCL_STR += "# Centers of mass\n"
TCL_STR += 'set ligninhead_com [measure center $ligninhead weight mass]\n'
TCL_STR += 'set lignintail_com [measure center $lignintail weight mass]\n'
TCL_STR += 'set bcdsecondary_com [measure center $bcdsecondary weight mass]\n'
TCL_STR += 'set bcd_com [measure center $bcd weight mass]\n\n'

TCL_STR += "# Vectors\n"
TCL_STR += "set head_tail_vector [vecsub $lignintail_com $ligninhead_com]\n"
TCL_STR += "set bcd_secondary_vector [vecsub $bcdsecondary_com $bcd_com]\n"
TCL_STR += "set head_tail_length [veclength $head_tail_vector]\n"
TCL_STR += "set bcd_secondary_length [veclength $bcd_secondary_vector]\n"
TCL_STR += "set head_tail_scale [expr $arrow_length / $head_tail_length]\n"
TCL_STR += "set bcd_secondary_scale [expr $arrow_length / $bcd_secondary_length]\n"
TCL_STR += "set ht_arrow_vector [vecscale $head_tail_scale $head_tail_vector]\n"
TCL_STR += "set bs_arrow_vector [vecscale $bcd_secondary_scale $bcd_secondary_vector]\n"
TCL_STR += "set ht_head_vector [vecscale 0.1 $ht_arrow_vector]\n"
TCL_STR += "set bs_head_vector [vecscale 0.1 $bs_arrow_vector]\n\n"

TCL_STR += "# Arrow & head positions \n"
TCL_STR += "set ht_arrow_end [vecadd $ligninhead_com $ht_arrow_vector]\n"
TCL_STR += "set ht_head_end [vecadd $ht_arrow_end $ht_head_vector]\n"
TCL_STR += "set bs_arrow_end [vecadd $bcd_com $bs_arrow_vector]\n"
TCL_STR += "set bs_head_end [vecadd $bs_arrow_end $bs_head_vector]\n\n"

TCL_STR += "# Draw arrows & heads\n"
TCL_STR += "draw color blue\n"
TCL_STR += "graphics top cylinder $ligninhead_com $ht_arrow_end radius 0.1 resolution 50\n"
TCL_STR += "graphics top cone $ht_arrow_end $ht_head_end radius 0.2 resolution 50\n"
TCL_STR += "graphics top cylinder $bcd_com $bs_arrow_end radius 0.1 resolution 50\n"
TCL_STR += "graphics top cone $bs_arrow_end $bs_head_end radius 0.2 resolution 50\n\n"

TCL_STR += "# Turn off axes\n"
TCL_STR += "axes location off\n\n"

TCL_STR += "set vector1 [vecsub $lignintail_com $ligninhead_com]\n"
TCL_STR += "set vector2 [vecsub $bcdsecondary_com $bcd_com]\n"
TCL_STR += "set dot_product [vecdot $vector1 $vector2]\n"
TCL_STR += 'puts "Dot product of the two vectors:"\n'
TCL_STR += "puts $dot_product\n"

with open(args.tcl, "w", encoding="utf-8") as f:
    f.write(TCL_STR)
