from multiprocessing import JoinableQueue, Process
import time
import os

import gromacs

# TODO how should the Analysis object communicate with the Universe?
class Analyzer(object):
	"""docstring for Analysis"""
	def __init__(self, filename, queue):
		# list of analysis objects
		self.__tasks = []
		self.__shared_queue = queue
	
	def run(self):
		# start a queue of size max 8, block if no empty slots
		for i in range(0, len(self.__tasks)):
			p = Process(target=__worker)
			p.start()
		
		for analysis in self.__queue:
			self.__shared_queue.put(analysis, True, None)
		self.__shared_queue.join()
	
	def add(self, analysis):
		"""docstring for attach_analysis"""
		self.__queue.append(analysis)
	
	def remove(self, analysis):
		"""docstring for remove_analysis"""
		self.__queue.append(analysis)

	def __worker(self):
		while True:
			item = q.get()
			item.run()
			q.task_done()
		
		
class Energy(Analysis):
	def __init__(self):
		super(Energy, self).__init__()

	def run(self):
		return self.__temperature()
	

class ContactMap(Analysis):
	def __init__(self, arg):
		super(ContactMap, self).__init__()
		self.arg = arg
		
	def run(self):
		pass

# TODO Implement
# class Hydration(Analysis):
# 	pass

# class Eed(Analysis):
# 	def run(self):
# 		# runs eed analysis
