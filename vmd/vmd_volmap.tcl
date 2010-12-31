


#load a starting gro file
#load in all the trajectories using a series of mol addfiles 
#do a rmsd alignment based on backbone
#output a single frame (last frame)
#compute the volume map 

#example
#mol new klv_nosol.gro
#mol addfile traj/traj.xtc
#volmap occupancy -allframes -combine avg -res 1.0 -o traj_rep.dx
