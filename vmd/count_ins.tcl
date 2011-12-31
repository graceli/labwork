# for a trajectory of inositol and some protein where the number of bound inositols is more than 8
# output the selection protein with the bound inositols (x angstroms away) as a pdb file

set numframes [molinfo top get numframes]

for {set i 0} {$i < $numframes} {incr i} {
	set around [atomselect top "same residue as resname INS and within 5 of protein" frame $i]
	set snap [atomselect top "not resname INS or (same residue as resname INS and within 5 of protein)" frame $i]
	set natoms [$around num]
	set nmols [expr $natoms/24]
	if { $nmols > 8 } {
		$snap writepdb output${i}_${nmols}.pdb
	}
}

close $outfile

