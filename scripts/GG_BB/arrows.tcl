# Arrow & head length
set arrow_length 9.0

# Atom selections
set ligninhead [atomselect top "index 171 173 175 181 184 186"]
set lignintail [atomselect top "index 147 149 151 157 160 162"]
set bcdsecondary [atomselect top "index 7 11 28 32 49 53 70 74 91 95 112 116 133 137"]
set bcd [atomselect top "resname BCDC"]

# Centers of mass
set ligninhead_com [measure center $ligninhead weight mass]
set lignintail_com [measure center $lignintail weight mass]
set bcdsecondary_com [measure center $bcdsecondary weight mass]
set bcd_com [measure center $bcd weight mass]

# Vectors
set head_tail_vector [vecsub $lignintail_com $ligninhead_com]
set bcd_secondary_vector [vecsub $bcdsecondary_com $bcd_com]
set head_tail_length [veclength $head_tail_vector]
set bcd_secondary_length [veclength $bcd_secondary_vector]
set head_tail_scale [expr $arrow_length / $head_tail_length]
set bcd_secondary_scale [expr $arrow_length / $bcd_secondary_length]
set ht_arrow_vector [vecscale $head_tail_scale $head_tail_vector]
set bs_arrow_vector [vecscale $bcd_secondary_scale $bcd_secondary_vector]
set ht_head_vector [vecscale 0.1 $ht_arrow_vector]
set bs_head_vector [vecscale 0.1 $bs_arrow_vector]

# Arrow & head positions 
set ht_arrow_end [vecadd $ligninhead_com $ht_arrow_vector]
set ht_head_end [vecadd $ht_arrow_end $ht_head_vector]
set bs_arrow_end [vecadd $bcd_com $bs_arrow_vector]
set bs_head_end [vecadd $bs_arrow_end $bs_head_vector]

# Draw arrows & heads
draw color blue
graphics top cylinder $ligninhead_com $ht_arrow_end radius 0.1 resolution 50
graphics top cone $ht_arrow_end $ht_head_end radius 0.2 resolution 50
graphics top cylinder $bcd_com $bs_arrow_end radius 0.1 resolution 50
graphics top cone $bs_arrow_end $bs_head_end radius 0.2 resolution 50

# Turn off axes
axes location off

set vector1 [vecsub $lignintail_com $ligninhead_com]
set vector2 [vecsub $bcdsecondary_com $bcd_com]
set dot_product [vecdot $vector1 $vector2]
puts "Dot product of the two vectors:"
puts $dot_product
