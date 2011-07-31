#!/usr/bin/python
import pickle
import os
import sys
import glob
import logging
import tarfile
import gromacs

class FileManager(object):
	"""
		Class to keep track of files inside a data directory
	"""
	def __init__(self, dir):
		print >> sys.stderr, dir, " is now tracked"
		# dictionary containing the processed state of the files in the directory
		self.__files = {}
		self.__persistent_storage = 'file_manager.pkl'
		self.__initialize_state(dir)
		
	def __del__(self):
		self.__pfile = open(self.__persistent_storage, 'w')
		pickle.dump(self.__files, self.__pfile)

	def processed_files(self):	
		return [ key for key in self.__files if self.__files[key] ]
		
	def unprocessed_files(self):
		return [ key for key in self.__files if not self.__files[key] ]
			
	def is_processed(self, filename):
		return self.__files[filename]
	
	def processed(self, filename):
		"""docstring for process"""
		self.__files[filename] = True
		
	def reset(self):
		"""docstring for reset"""
		for key in unprocessed:
			self.__files[key] = False				
	
	def __initialize_state(self, dir):
		""" record the initial state of the directory 
			at first time the object is initialized 
			files added to the directory will not be seen
		"""
		if os.path.exists(self.__persistent_storage):
			self.__pfile = open(self.__persistent_storage, 'r')
			self.__files = pickle.load(self.__pfile)
			self.__pfile.close()
		else:
			self.__pfile = open(self.__persistent_storage, 'w')

		if len(self.__files) == 0:
			# no previous persistence, start fresh
			files_list = os.listdir(dir)
			for f in files_list:
				self.__files[f] = False
			pickle.dump(self.__files, self.__pfile)
			self.__pfile.close()

class SH3FileSystem(object):
	"""docstring for SH3Filesystem"""
	def __init__(self, root):
		self.root = root
		self.fm = FileManager(root)
	
	def xtc_files(self):
		for tar in self.fm.unprocessed_files():
			tar_abs = os.path.join(self.root, tar)
			print tar_abs
			tar_manager = SH3Tarfile(tar_abs, '/dev/shm')
			self.__index(tar_manager)
			processed_files = tar_manager.index()
			yield processed_files
			
			# Mark tarfile processed if all xtcs extracted were indexed
			if tar_manager.done():
				tar_manager.processed(tar_abs)

	def __index(self, tar_manager):
		tar_manager.index()

	def num_tarfiles(self):
		return len(self.files)
	
class SH3Tarfile():
	def __init__(self, filename, source, tempdest='/dev/shm', dest='results.h5'):
		logging.basicConfig(filename='sh3tarfile.log',level=logging.DEBUG)
		self.__tarfile = filename
		self.__location = source
		self.__output = dest
		self.__scratch = tempdest
		self.__num_index = 0
		self.__num_traj = 0
	
	def index(self):
		if not tarfile.is_tarfile(self.__tarfile):
			return None
		
		# open tarfile with no compression
		self.__unpack_tarfile()
		trajectories = glob.glob("/dev/shm/*.xtc")
		self.__num_traj = len(trajectories)
		logging.info("%d trajectories found", self.__num_traj)
		
		processed = []
		for xtc in trajectories:
			print "xtc", xtc
			#parse replica and sequence number
			basename, replica_num, sequence = self.__parse_name(xtc)
			print "basename", basename

			# calculate average temperature for the small xtc
			temp = self.__temperature(self.__temp_path(basename))

			# TODO: write a row to the h5 file
			print f, basename, replica_num, sequence, temp
			self.__num_indexed += 1
			processed.append(xtc)
		
		return processed

	def done(self):
		# returns the number of files in the tar file
		return self.__num_indexed == self.__num_traj
	
	def num_xtc(self):
		"""docstring for num_xtc"""
		return self.__num_traj
		
	def num_edr(self):
		return self.__num_edr

	def __temp_path(filename):
		return os.path.join(self.__scratch, filename)
	
	def __source_path(filename):
		return os.path.join(self.__location, filename)
	
	def __output_path(filename):
		return os.path.join(self.__output, filename)

	def __unpack_tarfile(self):
		tf_handle = tarfile.open(os.path.join(self.__location, self.__tarfile), 'r:')
		tf_handle.extractall("/dev/shm")

	
	def __get_base_xtc_name(self, name):
		filename = os.path.basename(name)
		basename,ext = os.path.splitext(filename)
		first, second = basename.split('_')
		return first
	
	def __parse_name(self, xtc):
		name = self.__get_base_xtc_name(xtc)
		rhalf, shalf = name.split('.')
		replica_num = rhalf[3:]
		sequence = int(shalf)
		return name,replica_num,sequence

	def __temperature(self, name):
		# find and process the edr file
		rc,stdout,stderr = gromacs.g_energy(input='Temperature', f=name + '.edr', stdout=False)
		temperature = float(stdout.split('\n')[-3].split()[1])
		return temperature
		
def main():
	"""docstring for main"""
	fs = SH3FileSystem('/project/pomes/grace/test/PRIOR_TO_RESTART_Wed_Oct_27_04:27:47_EDT_2010/output/data')
	for batch in fs.xtc_files():
		print "extracted batch"
		print batch
		sys.exit(0)
	
if __name__ == '__main__':
	main()
