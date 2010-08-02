import os
import subprocess
import shlex

class Xtc(object):
	"""docstring for XTC"""
	def __init__(self, root, dataroot, xtcfile, tprfile):		
		print xtcfile

		prefix = 'ST'
		self.root = root	
		
		#make sure the path parts are chopped off
		filename = os.path.basename(xtcfile)
		basename = filename[0:len(filename)-4]

		self.basename = basename

		noprebasename = basename[len(prefix):]
		print noprebasename

		self.replicanum, self.seqnum, self.temp = noprebasename.split('.')
		self.dataroot = dataroot
		
		edr = basename + '.edr'
		xtc = basename + '.xtc'
				
		self.tprname = os.path.join(self.root, tprfile)
		self.xtcname = os.path.join(os.path.join(dataroot, 'xtcs'), xtc)
		self.edrname = os.path.join(os.path.join(dataroot, 'edr'), edr) 

		print "initialized for analysis"
		print self.tprname, self.xtcname, self.edrname, self.dataroot

	def meta():
		return {'replicanum':self.replicanum, 'seqnum':self.seqnum,'temp':self.temp}	
	
	def _createOutput(self, analysisName):
		output = os.path.join(self.root, analysisName)
		if not os.path.exists(output):
			os.mkdir(os.path.join(self.root, analysisName))
		return output

	def rg(self):
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


	def eed(self):
		eedpath = self._createOutput('eed')
		indexfile = os.path.join(self.root, 'eed.ndx')
		selection = {'group1' : 0, 'group2' : 1}
		group1 = selection['group1']
		group2 = selection['group2']
		command = "g_eed -f %s -s %s -n %s -o %s" % (self.xtcname,self.tprname,indexfile,os.path.join(eedpath,self.basename))
		code = subprocess.Popen(shlex.split(command),  stdin=subprocess.PIPE).communicate("%s %s" % (selection['group1'], selection['group2']))	

		return eedpath	
	
	def rama(self):
		ramapath = self._createOutput('rama')
		indexfile = os.path.join(self.root, 'rama.ndx')
		command = "g_rama -f %s -s %s -o %s" % (self.xtcname,self.tprname,os.path.join(ramapath,self.basename))
		code = subprocess.Popen(shlex.split(command), stdout=open(os.devnull, 'w'), stderr=open(os.devnull,'w'))
		return ramapath  		
		

	# add an analysis call here	
	#def dihedral(**selection):
		#mkdir dihedral
		#class g_rama ..
		#indexfile = 'dihedral.ndx'
		#selection = {'group1':'Protein'}
	#	return path
