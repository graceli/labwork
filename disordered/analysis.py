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
				self.__loader.load(analysis)
				self.__task_queue.task_done()

class Analysis(object):
	def __init__(self, location='/dev/shm', name='analysis', table_type=None, cols=0):
		self.__location = location
		self.__analysis_name = name
		self.__table_type = table_type
		self.__num_columns = cols
	
	def run(self, xtc='', tpr=''):
		pass
	
	def name(self):
		return self.__analysis_name
	
	def files(self):
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
	
class ContactMap(Analysis):
	def __init__(self):
		super(ContactMap, self).__init__(name='contact_map', table_type=rowtypes.ContactMapTable, cols=59)

	def files(self, xtc):
		return [ xtc + '.dist.txt', xtc + '.contact.txt', xtc + '.q.txt' ]

	def types(self):
		return [ rowtypes.ContactMapTable, rowtypes.ContactMapTable, rowtypes.QTable ]

	def type_names(self):
		return ['contact_map_table', 'contact_map_table', 'q_table']

	def run(self, xtc='', tpr='sh3.tpr'):
		super(ContactMap, self).run(xtc, tpr)
		self.__extract(xtc, tpr)
		base,ext = os.path.splitext(xtc)
		
	def __extract(self, xtc, tpr):		
		name = os.path.join(self.location(), xtc)
	
		print "PID", os.getpid(), "is extracting:", xtc, tpr, "to", name

		selection = {'group1':'C-alpha'}
		command = "/home/grace/bin/g_mdmat_g -f %s -s %s -t 0.6 -mean %s -txt-dist %s.dist.txt -txt-contact %s.contact.txt -txt-native %s.q.txt" % (xtc, tpr, name, name, name, name)
		(stdout, stderr) = subprocess.Popen(shlex.split(command),  stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE).communicate("%s" %(selection['group1']))


def main():
	print "testing ..."
	contact_test = ContactMap()
	print "files", contact_test.files('test.xtc')
	print "types", contact_test.type_names()
	contact_test.run(xtc='test.xtc', tpr='sh3_native.tpr')

if __name__ == '__main__':
	main()
