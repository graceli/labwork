import subprocess
import shlex
from optparse import OptionParser

usage = """ %prog <protein gro> <solvent gro> <topology> """
parser = OptionParser(usage)
parser.add_option("-o", dest="outputFile", help="name of the output file")
options, args = parser.parse_args()


if len(args) < 3:
	parser.error("not enough input files")
	sys.exit(0)

protein = args[0]
solvent = args[1]
top = args[2]

if not os.path.exists(protein):
	print "Error: protein gro file was not found"
	sys.exit(0)
	
if not os.path.exists(solvent):
	print "Error: solvent file was not found"
	sys.exit(0)

if not os.path.exists(top):
	print "Error: topology file was not found"
	sys.exit(0)


command="genbox -cp %s -cs %s -p %s -o" % (protein, solvent, top, parser.outputFile)
code = subprocess.Popen(shlex(command).split)

if code:
	print "Info: Finished with code", code
