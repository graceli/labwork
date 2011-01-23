# attempt 0 at constructing a more reusable and more readable 
# post-processing script

import os
import glob
import sys
import logging

#if len(sys.argv) < 2:
#	print "please enter a path or an expression for a path"
#	sys.exit(1)

#globPath = sys.argv[1] + '/' + 'sys*.xtc';
#print "will glob path", globPath

#dirList = glob.glob(sys.argv[1]+ '/'+'sys*.xtc')

appName = sys.argv[0]
dirList = glob.glob("*/*sys*.xtc")
CENTERGRP = 15
OUTPUTGRP = 16 
TPRFN = 'sys0_2.tpr'
TPRNOSOL = 'nosol.tpr'

numFilesToProcess = len(dirList)
#configure logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='test.log',
                    filemode='w')

logging.info("processing %d xtc files", numFilesToProcess)

for fname in dirList:
	logging.info("converting %s ...", fname)
	basefn = os.path.splitext(os.path.basename(fname))[0]
	output_middle = basefn  + "_cnosol"
	
	#ideally pass parameters to the gromacs executable with out remembering the gromacs syntax or by remembering 
	#a simplified one, and the input gets logged 
	os.system('echo %(CENTERGRP)d %(OUTPUTGRP)d | trjconv -f %(fname)s -s %(TPRFN)s -pbc mol -center -n center.ndx -o %(output_middle)s' % vars())
	

centeredFilesList = glob.glob("*/*cnosol.xtc")
numFilesProcessed = len(centeredFilesList)
logging.info("%d files were centered and stripped of solvent", numFilesProcessed)

if numFilesToProcess == numFilesProcessed:
	logging.info("concatenating the converted xtcs into one xtc: %d files", numFilesProcessed)
	os.system("trjcat -f */*cnosol.xtc -cat -dt 2 -o GA4_beta_scyllo_dt2_for_volmap.xtc")

logging.error("did not finish converting all trajectories something went wrong.")

# ideally parse out the error messages and store somewhere so I know exactly what went on
# os.system('echo 1 0 | trjconv -f %(output_middle)s -s %(TPRNOSOL)s -fit rot+trans -o %(output_final)s' % vars())	
