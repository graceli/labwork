from MDAnalysis import *
u = Universe('out.gro', 'tester.xtc')
#print u.atoms
#u.coordinates()

test = u.selectAtoms("protein")
#test.radiusOfGyration()

for ts in u.trajectory:
	# note that test is linked to ts object that gets updated 
	# when ts object is updated
	print ts.frame, test.radiusOfGyration()
