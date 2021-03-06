# computes the volume occupancy map for a given trajectory with cosolvent and protein
set filename [lindex $argv 0]
set groname [lindex $argv 1]

set file [file rootname $filename]
mol load gro $groname xtc $filename
volmap occupancy [atomselect top "resname 0YB"] -allframes -combine avg -res 1.0 -o /dev/shm/grace/${file}_glcnac.dx
volmap occupancy [atomselect top "protein"] -allframes -combine avg -res 1.0 -o /dev/shm/grace/${file}_protein.dx
mol delete top
exit

