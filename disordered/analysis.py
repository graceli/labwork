from multiprocessing import JoinableQueue, Process
from Queue import Empty
import time
import os
import subprocess
import shlex
import sys

#import gromacs
import loader
import file_system

class Analyzer(object):
	def __init__(self, data_root, working_dir, tpr, index=True, index_output='index.h5'):
		# list of analysis objects
		self.__analyses = []
		self.__working_dir = working_dir
		self.__fs = file_system.SH3FileSystem(data_root, index=True, index_output=index_output)
		self.__loader = loader.Loader(working_dir)
		self.__task_queue = JoinableQueue(8)
		self.__tpr = tpr

	def run(self):
		# start a queue of size max 8, block if no empty slots
		
		# populate the task queue with (analysis, xtc) items 
		for batch in self.__fs.xtc_files():
			for xtc in batch:
				for analysis in self.__analyses:
					print "queuing", analysis.name(), "and", xtc
					self.__task_queue.put([analysis, xtc], True, None)

			print "waiting for these tasks to finish"
			self.__task_queue.join()
	
		for i in range(0, 8):
			p = Process(target=self.__worker)
			p.start()

	def add(self, analysis):
		self.__analyses.append(analysis)
	
	def remove(self, analysis):
		self.__analyses.append(analysis)

	def __worker(self):
		# get one item from queue
		# block if queue is empty
		print "worker process started"
		while True:
			try:
				# timeout after 30 seconds
				analysis,xtc = self.__task_queue.get(True, 30)
			except Empty:
				break
			else:
				analysis.run(xtc, self.__tpr)
				self.__task_queue.task_done()

	# def load(self):
	# 	self.__aloader.load(base + '.xvg', 'rg', rowtypes.RGTable, 3)
	# 	self.__aloader.load(base + '.xvg', 'sas', rowtypes.SASTable, 3)
	# 	self.__aloader.load(base + '.xvg', 'eed', rowtypes.EETable, 3)
	# 	self.__aloader.load(base + '.xvg', 'rama', rowtypes.RamaTable, 3)
	# 	self.__aloader.load(base + '.q.txt', 'mdmat.q', rowtypes.QTable, 3)
	# 	os.system("rm \#*")

class Analysis(object):
	output = 'analysis.h5'
	
	def __init__(self, location='/dev/shm', name='analysis', table_type=None, cols=0):
		self.__location = location
		self.__analysis_name = name
		self.__table_type = table_type
		self.__num_columns = cols
		self.__files = []
	
	def run(self, xtc='', tpr=''):
		pass
	
	def name(self):
		return self.__analysis_name
	
	def files(self):
		return self.__files
	
	def num_columns(self):
		return self.__num_columns

	def table_structure(self):
		return self.__table_type
	
	def _make_output_path(self, analysisName):
		output = os.path.join(self.__location, self.__analysis_name)
		if not os.path.exists(output):
			os.mkdir(output)
		return output
	
class ContactMap(Analysis):
	def __init__(self):
		super(ContactMap, self).__init__(name='contact_map', table_type=rowtypes.ContactMapTable)
		
	def run(self, xtc='', tpr='sh3.tpr'):
		super(ContactMap, self).run(xtc, tpr)
		# extract
		self.__extract(xtc, tpr)
		base,ext = os.path.splitext(xtc)
		#load
		#self.__aloader.load(base + '.contact.txt', 'mdmat.contact', rowtypes.ContactMapTable, 3)
		
	def __extract(self, xtc, tpr):
		path = self._make_output_path('mdmat')
		name = os.path.join(path, xtc)
	
		print "PID", os.getpid(), "is extracting:", xtc, tpr, "to", name

		selection = {'group1':'C-alpha'}

		command = "/home/grace/bin/g_mdmat_g -f %s -s %s -t 0.6 -mean %s -txt-dist %s.dist.txt -txt-contact %s.contact.txt -txt-native %s.q.txt" % (xtc, tpr, name, name, name, name)
		(stdout, stderr) = subprocess.Popen(shlex.split(command),  stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE).communicate("%s" %(selection['group1']))

		#print stderr		
		self.__files = [ name + '.txt', name + '.dist.txt', name + '.contact.txt', name + '.q.txt' ]

		return path





