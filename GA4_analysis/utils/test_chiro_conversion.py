# attempt 0 at constructing a more reusable and more readable 
# post-processing script


#TODO: factor out the required files to run the script

import os
import glob
import sys
import logging
import getopt
import fnmatch

def setupLogger(appName):
	logger = logging.getLogger(appName)
	logger.setLevel(logging.DEBUG)
	
	ch = logging.StreamHandler()
	fh = logging.FileHandler(appName + '.log')

	#set output levels for each handler
	ch.setLevel(logging.DEBUG)
	fh.setLevel(logging.DEBUG)
	
	formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
	#set formatting for the handlers
	ch.setFormatter(formatter)
	fh.setFormatter(formatter)

	#attach handlers to logger
	logger.addHandler(ch)
	logger.addHandler(fh)
	
	return logger


def locate(pattern, root=os.curdir):
	filesList = []
	for path, dirs, files in os.walk(os.path.abspath(root)):
		#filesList = filesList + os.path.join(path, fnmatch.filter(files, pattern))
			# note that yield is just a iterator
			# while called again in the loop, it picks up where it last left off
			#yield os.path.join(path, filename)
			# but in my case, I just want this function to create a giant list
		for filename in fnmatch.filter(files,pattern):
			filesList.append(os.path.join(path, filename))

	return filesList
	
appName = sys.argv[0]
logger = setupLogger(appName)
dirList = locate("*.xtc") #glob.glob("*/*sys*.xtc")

#print dirList, len(dirList)


CENTERGRP = 15
OUTPUTGRP = 16 
TPRFN = 'sys0_2.tpr'
TPRNOSOL = 'nosol.tpr'

numFilesToProcess = len(dirList)
#configure logging
#logging.basicConfig(level=logging.DEBUG,
#                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
#                    datefmt='%m-%d %H:%M',
#                    filename='test.log',
#                    filemode='w')

logger.info("processing %d xtc files", numFilesToProcess)

if numFilesToProcess == 0:
	logger.info("No files to process ... quitting")
	sys.exit(0)

for fname in dirList:
	logger.info("converting %s ...", fname)
	basefn = os.path.splitext(os.path.basename(fname))[0]
	output_middle = basefn  + "_cnosol"
	
	#ideally pass parameters to the gromacs executable with out remembering the gromacs syntax or by remembering 
	#a simplified one, and the input gets logged 
	os.system('echo %(CENTERGRP)d %(OUTPUTGRP)d | trjconv -f %(fname)s -s %(TPRFN)s -pbc mol -center -n center.ndx -o %(output_middle)s' % vars())
	

centeredFilesList = glob.glob("*/*cnosol.xtc")
numFilesProcessed = len(centeredFilesList)
logger.info("%d files were centered and stripped of solvent", numFilesProcessed)

if numFilesToProcess == numFilesProcessed:
	logger.info("concatenating the converted xtcs into one xtc: %d files", numFilesProcessed)
	os.system("trjcat -f */*cnosol.xtc -cat -dt 2 -o GA4_beta_chiro_dt2_for_volmap.xtc")
else:
	logger.error("did not finish converting all trajectories something went wrong.")

# ideally parse out the error messages and store somewhere so I know exactly what went on
os.system('echo 1 0 | trjconv -f GA4_beta_chiro_dt2_for_volmap.xtc -s %(TPRNOSOL)s -fit rot+trans -o GA4_beta_chiro_dt2_for_volmap_fit.xtc' % vars())	
