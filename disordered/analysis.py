from multiprocessing import JoinableQueue, Process
import time
import os
#import gromacs
import loader
import file_system

class Analyzer(object):
	def __init__(self, data_root, working_dir, tpr, index=True, index_output='index.h5'):
		# list of analysis objects
		self.__analyses = []
		self.__working_dir = working_dir
		self.__fs = file_system.SH3FileSystem(data_root, index=True, index_output=index_output)
		self.__task_queue = JoinableQueue(8)
		self.__tpr = tpr

	def run(self):
		# start a queue of size max 8, block if no empty slots
		for i in range(0, 8):
			p = Process(target=self.__worker)
			p.start()
		
		# populate the task queue with (analysis, xtc) items 
		for batch in self.__fs.xtc_files():
			for analysis in self.__analyses:
				for a,xtc in zip(analysis, batch):
					print "queuing", a.name(), "and", xtc
					self.__task_queue.put([a, xtc], True, None)

			print "waiting for these tasks to finish"
			self.__shared_queue.join()
	
	def add(self, analysis):
		self.__analyses.append(analysis)
	
	def remove(self, analysis):
		self.__analyses.append(analysis)

	def __worker(self):
		# get one item from queue
		# block if queue is empty
		print "worker process started"
		while True:
			analysis,xtc = self.__task_queue.get()
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
	
	def __init__(self, location='/dev/shm', name='analysis'):
		self.__location = location
		self.__analysis_name = name
		self.__aloader = loader.Loader(location, location)
	
	def run(self):
		pass
	
	def _make_output_path(self, analysisName):
		output = os.path.join(self.__location, self.__analysis_name)
		if not os.path.exists(output):
			os.mkdir(output)
		return output
	
class ContactMap(Analysis):
	def __init__(self):
		super(ContactMap, self).__init__(name='contact_map')
		
	def run(self, xtc):
		super(ContactMap, self).run()
		# extract
		self.__extract(xtc)
		print xtc
		base,ext = os.path.splitext(xtc)
		#load
		self.__aloader.load(base + '.contact.txt', 'mdmat.contact', rowtypes.ContactMapTable, 3)
		
	def __extract(self, xtc):
		path = self._make_output_path('mdmat')
		indexfile = os.path.join(self.root, 'mdmat.ndx')
		name = os.path.join(mdmatpath,self.basename)
		selection = {'group1':'C-alpha'}
		command = "/home/grace/bin/g_mdmat_g -f %s -s %s -t 0.6 -mean %s -txt-dist %s.dist.txt -txt-contact %s.contact.txt -txt-native %s.q.txt" % (xtc, self.tprname, name, name, name, name)
		code = subprocess.Popen(shlex.split(command),  stdout=open(os.devnull, 'w'), stderr=open(os.devnull,'w'), stdin=subprocess.PIPE).communicate("%s" %(selection['group1']))
		return path





