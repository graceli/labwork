
set selectionStrAssembly "protein and not hydrogen or resname INS"
set selectionProtein "protein and not hydrogen"
#set numRes 186

#timestep in the trajectory - so that data is reported vs time rather than frame
set timestep 1

set outfile [open test.txt w]
set numframes [molinfo top get numframes]
puts $numframes

for {set i 0} {$i < $numframes} {incr i} {
	set currtime [expr $i * $timestep ]
	puts -nonewline $outfile "$currtime "
	

	set myAssembly [atomselect top $selectionProtein frame $i]
                set SASComponent [measure sasa 1.4 $myAssembly]
                puts -nonewline $outfile "$SASComponent "


                #hydrophobic subassembly
                set selectionStrComponent "sidechain and not hydrogen and not oxygen and not nitrogen"
                #puts "$currtime: residue $j"
                set myComponent [atomselect top $selectionStrComponent frame $i]
                #calculate surface area for current residue
                set SASComponent [measure sasa 1.4 $myAssembly -restrict $myComponent]
                puts -nonewline $outfile "$SASComponent "

                #hydrophobic subassembly
                set selectionStrComponent "(resname LYSH and name NZ) or (resname GLU and name CD OE1 OE2)"
                #puts "$currtime: residue $j"
                set myComponent [atomselect top $selectionStrComponent frame $i]
                #calculate surface area for current residue
                set SASComponent [measure sasa 1.4 $myAssembly -restrict $myComponent]
                puts -nonewline $outfile "$SASComponent "
	
	

	set myAssembly [atomselect top $selectionStrAssembly frame $i]
		set SASComponent [measure sasa 1.4 $myAssembly]
                puts -nonewline $outfile "$SASComponent "
	
	
		#hydrophobic subassembly
		set selectionStrComponent "sidechain and not hydrogen and not oxygen and not nitrogen"
		#puts "$currtime: residue $j"
		set myComponent [atomselect top $selectionStrComponent frame $i]
		#calculate surface area for current residue
		set SASComponent [measure sasa 1.4 $myAssembly -restrict $myComponent]
		puts -nonewline $outfile "$SASComponent "

		#hydrophobic subassembly
		set selectionStrComponent "(resname LYSH and name NZ) or (resname GLU and name CD OE1 OE2)"
		#puts "$currtime: residue $j"
		set myComponent [atomselect top $selectionStrComponent frame $i]
		#calculate surface area for current residue
		set SASComponent [measure sasa 1.4 $myAssembly -restrict $myComponent]
		puts -nonewline $outfile "$SASComponent "
	puts $outfile " ";
}
close $outfile

exit
