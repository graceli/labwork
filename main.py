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

	disklocation = os.getcwd()	
	
	if(USE_DEVSHM == True):
		templocation = '/dev/shm'
	else:
		templocation = disklocation

	aloader = Loader(templocation)

	#list tar filse
	tarfileslist = glob.glob("*STDR_running*.tar") 
	assert len(tarfileslist) > 0, "there are no tar files in this directory"
	
	print "these are the tar files to be analyzed", tarfileslist

	fnull = open(os.devnull, 'w')

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

		# look for xtc/ dir and edr dir
		#command = "ls -l %(templocation)s/STDR_running/xtcs %(templocation)s/STDR_running/edr" % vars()
		#subprocess.check_call(shlex.split(command), stdout=fnull, stderr=fnull)
		
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
		print "analysis file is here", aloader._result.location
		os.system("cd /dev/shm; rm -rf STDR_running*; mkdir test; mv * test; tar cvf test.tar test/ --remove-files; cp test.tar %s; rm -rf /dev/shm/*" % disklocation)


if __name__ == '__main__':
	main()

