import os
import subprocess
import shlex

class Xtc(object):
	"""docstring for XTC"""
	def __init__(self, root, path, xtcfile, tprfile):		
		prefix = 'ST'
		self.root = root	
		basename = xtcfile[0:len(xtcfile)-4]
		self.basename = basename
		noprebasename = xtcfile[len(prefix):len(xtcfile)-4]
		self.replicanum, self.seqnum, self.temp = noprebasename.split('.')
		self.path = path
		
		edr = basename + '.edr'
		xtc = basename + '.xtc'
				
		self.tprname = os.path.join(self.root, tprfile)
		self.xtcname = os.path.join(os.path.join(path, 'xtcs'), xtc)
		self.edrname = os.path.join(os.path.join(path, 'edr'), edr) 

		print "initialized for analysis"
		print self.tprname, self.xtcname, self.edrname, self.path		
		#mkdir basename

	def meta():
		return {'replicanum':self.replicanum, 'seqnum':self.seqnum,'temp':self.temp}	
	
	def _createOutput(self, analysisName):
		output = os.path.join(self.root, analysisName)
		if not os.path.exists(output):
			os.mkdir(os.path.join(self.root, analysisName))
		return output

	def rg(self):
		#current dir/
		#    rg.ndx
		#    sasa.ndx
		#    xtc/
		#    edr/
		#    system.tpr
		
		output = self._createOutput('rg')
		rgpath = output	
		indexfile = os.path.join(self.root, 'rg.ndx')
		selection = {'group1' : 'Protein'}
		
		group = selection['group1']
		command = "g_gyrate -f %s -s %s -n %s -o %s" % (self.xtcname,self.tprname,indexfile,os.path.join(rgpath,self.basename))

		print "outputting in", output
		print output, "created"
		print "checking in", self.root
		print "checking for", indexfile
		print "checking for", selection['group1']
		print "running ", command

		code = subprocess.Popen(shlex.split(command), stdin=subprocess.PIPE).communicate(selection['group1'])

		return rgpath
		
	def sasa(self):
		output = self._createOutput('sas')
		saspath = output
		indexfile = os.path.join(self.root, 'sas.ndx')
		selection = {'group1':'Protein', 'group2':'Protein'}

		group1 = selection['group1']
		group2 = selection['group2']
		command = "g_sas -f %s -s %s -n %s -o %s" % (self.xtcname,self.tprname,indexfile,os.path.join(saspath,self.basename))
		code = subprocess.Popen(shlex.split(command),  stdin=subprocess.PIPE).communicate("%s %s" % (selection['group1'], selection['group2']))	

		return saspath
		
	#def dihedral(**selection):
		#mkdir dihedral
		#class g_rama ..
		#indexfile = 'dihedral.ndx'
		#selection = {'group1':'Protein'}
	#	return path
