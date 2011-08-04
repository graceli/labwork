from optparse import OptionParser
import os
import file_system 
import analysis

def working_path(filename):
	return os.path.join(os.getcwd(), filename)

def main():
	"""docstring for main"""
	# parser = OptionParser()
	# parser.add_option("-f", "--tarfile-root-dir", action="store_true", dest="root", default=False,
	# help="sets the root directory containing all the tarfiles with trajectories (xtcs) to analyze")
	# 															
	# (options, args) = parser.parse_args()
	# if len(args) < 1:
	# 	parser.error("Please specify a .h5 input file")

	data_root = '/project/pomes/grace/test/PRIOR_TO_RESTART_Wed_Oct_27_04:27:47_EDT_2010/output/data'

	# initialize analysis
	analyzer = analysis.Analyzer(data_root, os.getcwd(), 'sh3_native.tpr')

	# queue up analysis tasks
	analyzer.add(analysis.ContactMap())
	analyzer.run()

if __name__ == '__main__':
	main()


# Setup:
# ======
# External shell script loops over the list of directories
# calls the main.py with name of directory
# Files produced:
# For each directory processed, the following is produced:
# index.h5 file consisting a listing of all the xtc files processed
# analysis.h5 consisting of all analysis attached to the system
# Cleans /dev/shm and copies all files to current working directory
