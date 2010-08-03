#!/usr/bin/env python

from loader import *
import xtc
import sys
import glob
import subprocess
import shlex
import logging

def interactive():
	answer= raw_input("would you liked to continue?[y/n]")
	if answer == "n":
		sys.exit(0)

def start_logging(session, loglevel):
	logname = 'a.log'
	logging.basicConfig(filename=session+'.log', level=loglevel,
				     	format='%(asctime)s %(levelname)s %(message)s')	
	logging.debug("=== starting %(session)s extraction ===" % vars())



def main():
	"""docstring for main"""
	
	option = sys.argv[1]
	sessionname = sys.argv[2]
	disklocation = sys.argv[3] 
	USE_DEVSHM = True
	
	start_logging(sessionname, logging.CRITICAL)

	#logging.debug("this is a test")

	if(USE_DEVSHM == True):
		templocation = '/dev/shm'
	else:
		templocation = disklocation

	aloader = Loader(templocation)

	#list tar filse
	tarfileslist = glob.glob(os.path.join(disklocation, '*STDR_running*.tar'))
	assert len(tarfileslist) > 0, "there are no tar files in this directory"
	
	#print "these are the tar files to be analyzed", tarfileslist

	fnull = open(os.devnull, 'w')
	
	numprocessed = 1
	for tarfile in tarfileslist:
		logging.debug("inflating", tarfile, "in /dev/shm")
	
		if option == "-i":
			interactive()

		logging.debug("copying tar file to /dev/shm")
	
		command = "cp %(tarfile)s %(templocation)s" % vars()
		code = subprocess.call(shlex.split(command))
		
		tarfile = os.path.join(templocation, tarfile)

		file, ext = os.path.splitext(tarfile)

		command = "tar xvf %(tarfile)s -C %(templocation)s" % vars() 
		args = shlex.split(command)
		code = subprocess.call(args)
		assert code == 0, "something bad happened during untarring. Stopping."
	
	
		tmpdataroot = os.path.join(templocation, 'STDR_running')
		
		logging.debug("the data root is %s", tmpdataroot)

		xtcpath = os.path.join(tmpdataroot, 'xtcs')
		xtcList = glob.glob(os.path.join(xtcpath,'*.xtc'))
		logging.info("these are the xtcs to be analyzed %s", xtcList)
		assert len(xtcList) > 0, "there are no xtcs found"
		

		if option == "-i":
			interactive()
			
		# analyze xtc files
		for xtcfile in xtcList:
			traj = xtc.Xtc(templocation, tmpdataroot, xtcfile, 'sh3.tpr')

			base, ext = os.path.splitext(os.path.basename(xtcfile))
			

			ramapath = traj.rama()
			rgpath = traj.rg()
			sasapath = traj.sasa()
			eedpath = traj.eed()
			mdmatpath = traj.mdmat()

			aloader.load(base + '.xvg', 'rg', rowtypes.RGTable, 3)
			aloader.load(base + '.xvg', 'sas', rowtypes.SASTable, 3)
			aloader.load(base + '.xvg', 'eed', rowtypes.EETable, 3)
			aloader.load(base + '.xvg', 'rama', rowtypes.RamaTable, 3)
			aloader.load(base + '.q.txt', 'mdmat.q', rowtypes.QTable, 3)
			aloader.load(base + '.contact.txt', 'mdmat.contact', rowtypes.ContactMapTable, 3)
			os.system("rm \#*")

			#aloader.load('eed', rowtypes.EEDTable, replicaMeta)
			#aloader.load('dihedral', rowtypes.DihedralTable, replicaMeta)
			#aloader.load('energy', rowtypes.EnergyTable, replicaMeta)

			if option == "-i":
				interactive()
	
		### need to clean up /dev/shm here after each tar is processed ###	

		# run a bash scriptlet
		# add additional analysis to scriptlet
		os.system("cd /dev/shm; mkdir preserve; mv sh3.tpr *.ndx preserve; rm -rf STDR_running*; mkdir STDR_running_analyzed_%(numprocessed)s; mv mdmat* rama* rg* sas* eed* STDR_running_analyzed_%(numprocessed)s; cp preserve/* .; tar cvf STDR_running_analyzed_%(numprocessed)s.tar STDR_running_analyzed_%(numprocessed)s; gzip STDR_running_analyzed_%(numprocessed)s.tar; cp STDR_running_analyzed_%(numprocessed)s.tar.gz %(disklocation)s; rm -rf STDR_running*; rm -r preserve; cd %(disklocation)s" % vars())

		numprocessed += 1

	#end of the tar loop
	os.system("cp /dev/shm/*.h5 %(disklocation)s" % vars())
	os.system("rm \#*")



if __name__ == '__main__':
	main()
	os.system("rm -rf /dev/shm/*")

