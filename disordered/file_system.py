#!/usr/bin/python
import pickle
import os
import sys
import glob
import logging
import tarfile

import subprocess
import shlex

#import gromacs

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
	def __init__(self, root, index=True, index_output='index.h5'):
		self.__perform_indexing = index
		self.root = root
		self.fm = FileManager(root)
	
	def xtc_files(self):
		processed_files = []
		for tar in self.fm.unprocessed_files():
			tar_abs = os.path.join(self.root, tar)

			print tar_abs
			print self.root

			tar_manager = SH3Tarfile(tar_abs, self.root)
			processed_files = tar_manager.index(index=self.__perform_indexing)

			yield processed_files
			
			# Mark tarfile processed if all xtcs extracted were indexed
			# TODO: should mark individual xtc files?
			if tar_manager.done():
				self.fm.processed(tar_abs)
		yield processed_files

	def num_tarfiles(self):
		return len(self.files)
	
class SH3Tarfile():
	def __init__(self, filename, source, tempdest='/dev/shm', dest='results.h5'):
		logging.basicConfig(filename='sh3tarfile.log',level=logging.DEBUG)
		logging.info("initializing SH3Tarfile with filename=%s in source=%s", filename, source) 
		self.__tarfile = filename
		self.__location = source
		self.__output = dest
		self.__tempdest = tempdest
		self.__num_indexed = 0
		self.__num_traj = 0
	
	def index(self, index=True):
		logging.info("SH3Tarfile: indexing ...")

		if not tarfile.is_tarfile(self.__tarfile):
			return None
		
		# open tarfile with no compression
		self.__unpack_tarfile()
		trajectories = glob.glob(self.__tempdest + "/" + "*.xtc")
		self.__num_traj = len(trajectories)
		logging.info("%d trajectories found in %s", self.__num_traj, self.__tempdest)
	
		assert self.__num_traj > 0, "SH3Tarfile: No trajectories were found"

		if index:
			f=open(os.path.join(self.__tempdest, 'index.txt'), 'a')

		processed = []
		for xtc in trajectories:
			#parse replica and sequence number
			basename, replica_num, sequence = self.__parse_name(xtc)
			# calculate average temperature for the small xtc
			temp = self.__temperature(self.__temp_path(basename))

			# TODO: write a row to the h5 file
			if index:
				print >> f, basename, replica_num, sequence, temp

			self.__num_indexed += 1
			processed.append(SH3XtcFile(os.path.basename(xtc), self.__tempdest, replica_num, sequence, temp))
		
		return processed

	def done(self):
		# returns the number of files in the tar file
		return self.__num_indexed == self.__num_traj
	
	def num_xtc(self):
		"""docstring for num_xtc"""
		return self.__num_traj
		
	def num_edr(self):
		return self.__num_edr

	def __temp_path(self, filename):
		return os.path.join(self.__tempdest, filename)
	
	def __source_path(self, filename):
		return os.path.join(self.__location, filename)
	
	def __output_path(self, filename):
		return os.path.join(self.__output, filename)

	def __unpack_tarfile(self):
		tf_handle = tarfile.open(os.path.join(self.__location, self.__tarfile), 'r:')
		tf_handle.extractall(self.__tempdest)

	
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
		#rc,stdout,stderr = gromacs.g_energy(input='Temperature', f=name + '.edr', o=name, stdout=False)
		# use subprocess to grab average temperature instead ... messy read only home dir on compute nodes
		selection = {'group1':'Temperature'}
		logging.info('getting temperature for %s', name)
		command = "g_energy -f %s -o %s" % (name, self.__temp_path(name))

		(stdout, stderr) = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=open(os.devnull,'w'), stdin=subprocess.PIPE).communicate(selection['group1'])
		temperature = float(stdout.split('\n')[-3].split()[1])
		return temperature

class File(object):
	def __init__(self, name, location, type='txt'):
		self._name = name
		self._location = location
		self._type = type
	
	def type(self):
		return self._type

	def name(self):
		return self._name

	def location(self):
		return self._location
	
	def params(self):
		pass

	def full_path(self):
		return os.path.join(self._location, self._name)

class SH3XtcFile(File):
	def __init__(self, name, location, replica, sequence, temperature):
		super(SH3XtcFile, self).__init__(name, location) 
		self.__replica = replica
		self.__sequence = sequence
		self.__temperature = temperature

	def params(self):
		return [self.__replica, self.__sequence, self.__temperature]

def main():
	fs = SH3FileSystem('/project/pomes/grace/test/PRIOR_TO_RESTART_Wed_Oct_27_04:27:47_EDT_2010/output/data')
	for batch in fs.xtc_files():
		logging.info("extracted batch")
		print batch
		sys.exit(0)
	
if __name__ == '__main__':
	main()
