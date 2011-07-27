#!/usr/bin/python
import pickle
import os
import sys

"""
	Class to keep track of files inside a data directory
"""

class FileManager(object):
	# dictionary containing the processed state of the files in the directory
	
	def __init__(self, dir):
		print >> sys.stderr, dir, " is now tracked" 
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
				
def main():
	"""docstring for main"""
	import time
	fm = FileManager('.')
	for f in fm.unprocessed_files():
		print >> sys.stderr, "processing",f
		fm.process(f)
		time.sleep(10)
	
if __name__ == '__main__':
	main()
