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
import rowtypes

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
		for i in range(0, 8):
			p = Process(target=self.__worker)
			p.start()

		for batch in self.__fs.xtc_files():
			print "batch", batch
			for xtc in batch:
				for analysis in self.__analyses:
					print "queuing", analysis.name(), "and", xtc.name()
					self.__task_queue.put([analysis, xtc], True, None)

			print "waiting for these tasks to finish"
			self.__task_queue.join()
			print "tasks have finished"

			print "PID", os.getpid(), "loading analysis"
			for xtc in batch:
				for a in self.__analyses:
					self.__loader.load(a, xtc)	

	def add(self, analysis):
		self.__analyses.append(analysis)
	
	def remove(self, analysis):
		self.__analyses.append(analysis)

	def __worker(self):
		# TODO: use pool because it looks like the processes sometimes don't die if it fails
		# get one item from queue
		# block if queue is empty
		while True:
			try:
				# timeout after 30 seconds
				analysis,xtc = self.__task_queue.get(True, 30)
			except Empty:
				break
			else:
				analysis.run(xtc, self.__tpr)
				self.__task_queue.task_done()

class Analysis(object):
	def __init__(self, location='/dev/shm', name='analysis', table_type=None, cols=0):
		self.__location = location
		self.__analysis_name = name
		self.__table_type = table_type
		self.__num_columns = cols
	
	def run(self, xtc=None, tpr=''):
		pass	

	def name(self):
		return self.__analysis_name
	
	def files(self, xtc):
		pass
	
	def types(self):
		pass
	
	def location(self):
		return self.__location
	
	def num_columns(self):
		return self.__num_columns

	def table_structure(self):
		return self.__table_type
	
	def _make_output_path(self, analysisName):
		output = os.path.join(self.__location, self.__analysis_name)
		if not os.path.exists(output):
			os.mkdir(output)
		return output

class Gyration(Analysis):
	
	def files(self, xtc):
		return [ xtc.name() + '.xvg' ]
	
	def types(self):
		return [ rowtypes.RGTable ]
	
	def type_names(self):
		return ['rg']
	
	def xtc_file(self):
		return self.__file
		
	
	def run(self, xtc=None, tpr='sh3.tpr'):
		self.__extract(xtc, tpr)
		self.__file = xtc

	def __extract(self, xtc, tpr):		
		name = os.path.join(self.location(), xtc.name())

		print "PID", os.getpid(), "is extracting:", xtc.name(), tpr, "to", name

		selection = {'group1':'Protein'}
		command = "g_gyrate -f %s -s %s " % (xtc.full_path(), tpr, name, name, name, name)
		(stdout, stderr) = subprocess.Popen(shlex.split(command),  stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE).communicate("%s" %(selection['group1']))
		print "PID", os.getpid(), "done extracting"
	
	
	
class ContactMap(Analysis):
	def __init__(self):
		super(ContactMap, self).__init__(name='contact_map', table_type=rowtypes.ContactMapTable, cols=59)
		self.__file = None

	def files(self, xtc):
		return [ xtc.name() + '.contact.txt', xtc.name() + '.q.txt' ]

	def types(self):
		return [ rowtypes.ContactMapTable, rowtypes.QTable ]

	def type_names(self):
		return ['contact_map', 'q_native']

	def xtc_file(self):
		return self.__file
	
	def run(self, xtc=None, tpr='sh3.tpr'):
		super(ContactMap, self).run(xtc, tpr)
		self.__extract(xtc, tpr)
		self.__file = xtc
		# base,ext = os.path.splitext(xtc)
		
	def __extract(self, xtc, tpr):		
		name = os.path.join(self.location(), xtc.name())
	
		print "PID", os.getpid(), "is extracting:", xtc.name(), tpr, "to", name

		selection = {'group1':'C-alpha'}
		command = "/home/grace/bin/g_mdmat_g -f %s -s %s -t 0.6 -mean %s -txt-dist %s.dist.txt -txt-contact %s.contact.txt -txt-native %s.q.txt" % (xtc.full_path(), tpr, name, name, name, name)
		(stdout, stderr) = subprocess.Popen(shlex.split(command),  stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE).communicate("%s" %(selection['group1']))
		print "PID", os.getpid(), "done extracting"

def main():
	print "testing ..."
	contact_test = ContactMap()
	print "files", contact_test.files('test.xtc')
	print "types", contact_test.type_names()
	contact_test.run(xtc='test.xtc', tpr='sh3_native.tpr')

if __name__ == '__main__':
	main()
