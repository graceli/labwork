import glob
import os
import logging
import subprocess
import sys

logging.basicConfig(filename='prepare.log', level=logging.DEBUG)
logger = logging.getLogger('prepare')
		
def setup_directories(nfiles):
	for i in range(nfiles):
		print "Creating run dir", i
		os.system("mkdir %(i)d" % vars())
		os.system("cp -rp *_%(i)d_*.top *_%(i)d_*renum.gro %(i)d/" % vars())
		# rename files
		os.system("mv %(i)d/*.gro  %(i)d/system_%(i)d_final.gro" % vars())

def clean_up(num):
	for i in range(num):
		print "Removing directories %(i)d" %vars()
		os.system("rm -rf %(i)d" % vars())

def test_compile(num, stage="equil"):
	""" compiles the run initialization files for the systems """
	logging.info("Test compiling the tpr files for stage=%s", stage)
	mdp_file = "%(stage)s.mdp" % vars()
	
	for i in range(num):
		gro = "system_%(i)d_final.gro" % vars()
		top = "system_%(i)d_final.top" % vars()
		try:	
			command = "grompp_mpi -f params/%(mdp_file)s -c %(i)d/%(gro)s -p %(i)d/%(top)s -o %(i)d/%(stage)s.tpr -po %(i)d/used.mdp" % vars()
			logging.debug("compiling system = (%s, %s, %s)", mdp_file, gro, top)
			logging.debug(command)
			subprocess.check_call(command, shell=True, stderr=subprocess.STDOUT, stdout=open('gromacs.log', 'a'))
		except subprocess.CalledProcessError:
			print "EXCEPTION"
			logging.debug("Error occurred: %s", sys.exc_info()[0])
	
def main():
	nfiles = 44
	print "preparing", nfiles, "systems"
	setup_directories(nfiles)
	test_compile(nfiles, stage='em')
	#test_compile(nfiles, stage="production")	

if __name__ == '__main__':
	main()
