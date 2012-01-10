# a script to wrap around HPSS recall, off load, listing, mechanisms

from optparse import OptionParser
import os
import sys
import logging

# provide command line options to recall files
# template scripts
# allows generation of the script only
# allows generation and submission
# keeps a log of your retrievals
# generate a bash script and submits it to archive queue


class Script:
	bash_script = " "



parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="write report to FILE", metavar="FILE")
parser.add_option("-l", "--log", dest="logname", help="output a log", metavar="LOG")
(options, args) = parser.parse_args()

print options
print args

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(options.logname)
logger.info("this is a log!")


