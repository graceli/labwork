set outfile [open contacts.xvg w]
set numframes [molinfo top get numframes]


for {set i 0} {$i < $numframes} {incr i} {

set mysel [atomselect top "resid 140 to 164 and within 3.5 of resid 6 to 36" frame $i]
set atomlist   [$mysel get index]
set numContacts [$mysel num]
#puts "list of atoms within 3 of protein:"
#puts "$atomlist"
#puts $numContacts

puts $outfile "$i $numContacts"
}
close $outfile
