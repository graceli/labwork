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
	tarfileslist = glob.glob("*STDR_running*.tar") #glob.glob(os.path.join(current, '*.tar'))
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
		#command = "ls -l %(templocation)s/STDR_running/xtcs %(templocation)s/STDR_running/edr" % vars()
		#subprocess.check_call(shlex.split(command), stdout=fnull, stderr=fnull)
		
		xtcpath = os.path.join(path, 'xtcs')
		xtcList = glob.glob(os.path.join(xtcpath,'*.xtc'))
		print "these are the xtcs to be analyzed", xtcList
		assert len(xtcList) > 0, "there are no xtcs found"
		

		if sys.argv[1] == "-i":
			interactive()
			
		# analyze xtc files
		for xtcfile in xtcList:
			traj = xtc.Xtc(templocation, path, xtcfile, 'sh3.tpr')

			base, ext = os.path.splitext(xtcfile)
			
			rgpath = traj.rg()
			aloader.load(base + '.xvg', 'rg', rowtypes.RGTable, 3)
	
			sasapath = traj.sasa()
			aloader.load(base + '.xvg', 'sas', rowtypes.SASTable, 3)

			#aloader.load('eed', rowtypes.EEDTable, replicaMeta)
			#aloader.load('dihedral', rowtypes.DihedralTable, replicaMeta)
			#aloader.load('energy', rowtypes.EnergyTable, replicaMeta)

			if sys.argv[1] == "-i":
				interactive()
		
		cleanup(tarfile, templocation, current)	

def cleanup(tarfile, temp, disklocation):
		print "cleaning up ..."
		command = "rm -rf *STDR_running*"
		os.system(command)

		command = "mkdir %(tarfile)s_analysis" % vars()
		subprocess.check_call(shlex.split(command))

		command = "cp %(temp)s/* %(temp)s/%(tarfile)s_analysis" % vars()
		os.system(command)

		command = "tar cvf %(tarfile)s_analysis.tar %(tarfile)s_analysis" % vars()
		subprocess.check_call(shlex.split(command))

		command = "rm -rf %(tarfile)s_analysis" % vars()
		subprocess.check_call(shlex.split(command))
		
		command = "cp %(tarfile)s_analysis.tar %(disklocation)" % vars()
		subprocess.check_call(shlex.split(command))

if __name__ == '__main__':
	main()

