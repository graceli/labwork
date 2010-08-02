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

	current = os.getcwd()	
	
	if(USE_DEVSHM == True):
		templocation = '/dev/shm'
	else:
		templocation = current

	aloader = Loader(current)
	
	#list tar filse
	tarfileslist = glob.glob("*.tar") #glob.glob(os.path.join(current, '*.tar'))
	assert len(tarfileslist) > 0, "there are no tar files in this directory"
	
	print "these are the tar files to be analyzed", tarfileslist

	#copy tar file into /dev/shm and attempt to inflate tar in /dev/shm
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
	
	
		path = os.path.join(templocation, 'STDR_running')
		
		print "now in", path

		# look for xtc/ dir and edr dir
		command = "ls -l STDR_running/xtcs STDR_running/edr"
		code = subprocess.call(shlex.split(command), stdout=fnull, stderr=fnull)
		assert code == 0, "xtc or edr not found here"
		
		xtcpath = os.path.join(path, 'xtcs')
		xtcList = os.listdir(xtcpath)
		print "these are the xtcs to be analyzed", xtcList
		

		if sys.argv[1] == "-i":
			interactive()
			
		# analyze xtc files
		for xtcfile in xtcList:
			traj = xtc.Xtc(templocation, path, xtcfile, 'sh3.tpr')

			base, ext = os.path.splitext(xtcfile)
			
			#rgpath = traj.rg()
			#aloader.load('rg', rowtypes.RGTable, 3)
	
			sasapath = traj.sasa()
			aloader.load(base+'.xvg', 'sas', rowtypes.SASTable, 3)

			#aloader.load('eed', rowtypes.EEDTable, replicaMeta)
			#aloader.load('dihedral', rowtypes.DihedralTable, replicaMeta)
			#aloader.load('energy', rowtypes.EnergyTable, replicaMeta)

			if sys.argv[1] == "-i":
				interactive()
			

if __name__ == '__main__':
	main()

