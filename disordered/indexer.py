import file_manager
import tarfile

# Example final output of the pytable file:
# lgwXXXX.YYYY xtc_name replica_num sequence_num average T
# ...
# ...

def xtc_files(members):
	for tarinfo in members:
		if os.path.splitext(tarinfo.name)[1] == '.xtc':
			yield tarinfo

# open the directory to index
f = open('list')
dirs = f.readlines()
for d in dirs:
	print "in", d.rstrip()
	
	fm = file_manager.FileManager(d)
	for f in fm.unprocessed_files():
		f_abs = os.path.join(d, f)
		
		print "processing", f_abs
		
		if tarfile.is_tarfile(f_abs):
			# open tarfile with no compression
			tf_handle = tarfile.open(f_abs, 'r:')
			# extract only the xtcs to /dev/shm
			# tf_handle.extractall("/dev/shm", members=xtc_files(tf_handle))
			tf_handle.extractall("/dev/shm")
			
		#untar files
		# trajectories = ...
		# for xtc in trajectories:
			#parse replica and sequence number
			# calculate average temperature for the small xtc
			# write a row to the h5 file
			
			
