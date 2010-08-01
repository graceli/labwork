class Xtc(object):
	"""docstring for XTC"""
	def __init__(self, path, xtcfile, tprfile):		
		prefix = 'ST'
		
		basename = xtcfile[0:len(filename)-4]
		noprebasename = xtcfile[len(prefix):len(filename)-4]
		self.replicanum, self.seqnum, self.temp = noprebasename.split(delimiter)
		self.path = path
		
		edr = basename + '.edr'
		xtc = basename + '.xtc'
				
		self.tprname = os.path.join(path, tprfile)
		self.xtcname = os.path.join(os.path.join(path, 'xtc'), xtc)
		self.edrname = os.path.join(os.path.join(path, 'edr'), edr) 
		
		#mkdir basename
		
	def rg(**selection):
		#current dir/
		#    rg.ndx
		#    sasa.ndx
		#    xtc/
		#    edr/
		#    system.tpr
		
		indexfile = os.path.join(self.path, 'rg.ndx')
		selection = {'group1':'Protein'}
		path = os.path.join(path, 'rg')
		
		#better to do with the subprocess module?
		command = "echo %(selection['group1']) | g_gyrate -f %(xtcname) -s %(tprname) -n %(indexfile) -o %(rgpath)/%(basename)_rg"
		os.system( command % vars())
				
		return path
		
	def sasa(**selection):
		indexfile = os.path.join(self.path, 'sasa.ndx')
		selection = {'group1':'Protein', 'group2':'Protein'}
		path = os.path.join(path, 'sas')
		
		os.system("echo %(selection['group1']) %(selection['group2']) | g_sas -f %(xtcname) -s %(tprname) -n %(indexfile) -o %(path)/%(basename)_sas")

		return path
		
	#def dihedral(**selection):
		#mkdir dihedral
		#class g_rama ..
		#indexfile = 'dihedral.ndx'
		#selection = {'group1':'Protein'}
	#	return path