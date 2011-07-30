import file_system

# Example final output of the pytable file:
# lgwXXXX.YYYY xtc_name replica_num sequence_num average T
# ...
# ...

# Architecture of this analysis package is beginning to be similar to MDAnalysis
# object that all the other classes knows about (how to implement a universe)
# only one in existence

class Universe(object):
	"""docstring for MDSystem"""
	def __init__(self, root):
		self.__fs = file_system.SH3FileSystem(root)
		self.__analyzer = analysis.Analyzer()
		self.__shared_queue = JoinableQueue(8)
		
	def analyze(self):
		# analyze all trajectories in the Universe
		for i in range(0, 8):
			p = Process(target=__worker)
		
		for trajs in fs.xtc_files():
			for xtc in trajs:
				self.__shared_queue.put(xtc)
		
	def __worker(self):
		traj = self.__shared_queue.get()
		self.__analyzer.run(traj)
		self.__shared_queue.task_done()
