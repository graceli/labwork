#!/usr/bin/env python

from loader import *
import xtc
import sys
import glob
import subprocess
import shlex

def interactive():
	answer= raw_input("would you liked to continue?[y/n]")
	if answer == "n":
		sys.exit(0)


def main():
	"""docstring for main"""

	USE_DEVSHM = True

	disklocation = sys.argv[2] 
	
	if(USE_DEVSHM == True):
		templocation = '/dev/shm'
	else:
		templocation = disklocation

	aloader = Loader(templocation)

	#list tar filse
	tarfileslist = glob.glob(os.path.join(disklocation, '*STDR_running*.tar'))
	assert len(tarfileslist) > 0, "there are no tar files in this directory"
	
	print "these are the tar files to be analyzed", tarfileslist

	fnull = open(os.devnull, 'w')
	
	numprocessed = 1
	for tarfile in tarfileslist:
		print "inflating", tarfile, "in /dev/shm"
	
		if sys.argv[1] == "-i":
			interactive()

		print "copying tar file to /dev/shm"
	
		command = "cp %(tarfile)s %(templocation)s" % vars()
		code = subprocess.call(shlex.split(command))
		
		tarfile = os.path.join(templocation, tarfile)

		file, ext = os.path.splitext(tarfile)

		command = "tar xvf %(tarfile)s -C %(templocation)s" % vars() 
		args = shlex.split(command)
		code = subprocess.call(args)
		assert code == 0, "something bad happened during untarring. Stopping."
	
	
		tmpdataroot = os.path.join(templocation, 'STDR_running')
		
		print "the data root is", tmpdataroot

		xtcpath = os.path.join(tmpdataroot, 'xtcs')
		xtcList = glob.glob(os.path.join(xtcpath,'*.xtc'))
		print "these are the xtcs to be analyzed", xtcList
		assert len(xtcList) > 0, "there are no xtcs found"
		

		if sys.argv[1] == "-i":
			interactive()
			
		# analyze xtc files
		for xtcfile in xtcList:
			traj = xtc.Xtc(templocation, tmpdataroot, xtcfile, 'sh3.tpr')

			base, ext = os.path.splitext(os.path.basename(xtcfile))
			
			rgpath = traj.rg()
			aloader.load(base + '.xvg', 'rg', rowtypes.RGTable, 3)
	
			sasapath = traj.sasa()
			aloader.load(base + '.xvg', 'sas', rowtypes.SASTable, 3)

			#aloader.load('eed', rowtypes.EEDTable, replicaMeta)
			#aloader.load('dihedral', rowtypes.DihedralTable, replicaMeta)
			#aloader.load('energy', rowtypes.EnergyTable, replicaMeta)

			if sys.argv[1] == "-i":
				interactive()
	
		### need to clean up /dev/shm here after each tar is processed ###	

		# run a bash scriptlet
		os.system("cd /dev/shm; mkdir preserve; mv sh3.tpr *.ndx preserve; rm -rf STDR_running*; mkdir STDR_running_analyzed_%(numprocessed); mv * STDR_running_analyzed_%(numprocessed); cp STDR_running_analyzed_i%(numprocessed)/preserve/* .; tar cvf STDR_running_analyzed_%(numprocessed).tar STDR_running_analyzed_%(numprocessed)/ --remove-files; gzip STDR_running_analyzed_%(numprocessed).tar; cp STDR_running_analyzed_%(numprocessed).tar.gz %(disklocation); rm -rf /dev/shm/STDR_running*; rm -r preserve; cd %(disklocation)" % vars())

		numprocessed += 1

	os.system("cp /dev/shm/*.h5 %(disklocation)")


if __name__ == '__main__':
	main()
	os.system("rm -rf /dev/shm/*")

