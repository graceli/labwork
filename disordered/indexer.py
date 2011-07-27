import os
import sys
import glob
import tarfile

import file_manager

# Example final output of the pytable file:
# lgwXXXX.YYYY xtc_name replica_num sequence_num average T
# ...
# ...

def parse_name(xtc):
	filename = os.path.basename(xtc)
	basename,ext = os.path.splitext(filename)
	first, second = basename.split('_')
	rhalf, shalf = first.split('.')
	replica_num = rhalf[3:]
	sequence = int(shalf)
	return basename,replica_num,sequence

def xtc_files(members):
	for tarinfo in members:
		if os.path.splitext(tarinfo.name)[1] == '.xtc':
			yield tarinfo

# open the directory to index
f = open('list')
dirs = f.readlines()
for d in dirs:
	d = d.rstrip()
	print >> sys.stderr, "in", d
	
	fm = file_manager.FileManager(d)
	for f in fm.unprocessed_files():
		f_abs = os.path.join(d, f)
		print >> sys.stderr, "processing", f
		if tarfile.is_tarfile(f_abs):
			# open tarfile with no compression
			tf_handle = tarfile.open(f_abs, 'r:')
			tf_handle.extractall("/dev/shm")
			
		trajectories = glob.glob("/dev/shm/*.xtc")
		print >> sys.stderr, len(trajectories)

		for xtc in trajectories:
			#parse replica and sequence number
			basename, replica_num, sequence = parse_name(xtc)
			# calculate average temperature for the small xtc
			# TODO
			# write a row to the h5 file
			print f, basename, replica_num, sequence
			os.remove(xtc)		
			print >> sys.stderr, "removed", xtc

		fm.processed(f)	
