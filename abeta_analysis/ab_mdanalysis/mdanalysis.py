from MDAnalysis import *

# this bit of code reads in a gro file and a corresponding trajectory
# selects the protein from the trajectory
# and extracts the radius of gyration of the protein

#snapshot = 'test_data/ab.gro'
snapshot = 'test_data/ab_prot_ins.tpr'
tester_traj = 'test_data/ab_tester.xtc'

u = Universe(snapshot, tester_traj)
#print u.atoms
#u.coordinates()

test = u.selectAtoms("protein")
#test.radiusOfGyration()

for ts in u.trajectory:
	# note that test is linked to ts object that gets updated 
	# when ts object is updated
	print ts.frame*1000, test.radiusOfGyration()/10
