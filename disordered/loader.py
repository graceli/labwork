# -*- coding: utf-8 -*-
import os
import tables
import xvgfile
import rowtypes

class Result:
	#wrappers around the pytable 
	#group -- type of analysis - structural, global etc
	#becomes a group in pytables
	#quantityName -- eg. Rg, DSSP, etc
	
	def __init__(self, targetpath, **initparams):
		#creates the pytable file and initializes a single group
		options = {
			'tabletitle' : 'default pytable store',
			'groupName' : 'Protein', 
			'groupDesc' : 'properties', 
			'filename' : 'analysis.h5'
		}
		options.update(initparams)
		
		self._filename = options['filename'] 
		self.targetpath = targetpath
		self.location = os.path.join(targetpath, self._filename)
		self._h5file = tables.openFile(self.location, mode = "a", title=options['tabletitle'])
		
		if not self._h5file.root.__contains__(options['groupName']):
			self._h5file.createGroup(self._h5file.root, options['groupName'], options['groupDesc'])

		#print "a h5 file", self._filename, "has been created in", self.location
		#print "and a group with name",options['groupName'],"has been added"

	def __del__(self):
		self._h5file.close()	
		
	def addToTable(self, data, **tableinitparams):
		options = {
				'groupName' : 'Protein',
				'tableName' : 'defaultname',
				'tableTitle' : 'defaulttitle',
				'tableStruct' : rowtypes.DefaultTable
		}
		options.update(tableinitparams)
		#print "group to be added to is", options['groupName']
	
		assert self._h5file.root.__contains__(options['groupName']) == True, "group is not found add the group first"
			
		group = self._h5file.root._f_getChild(options['groupName'])

		if group.__contains__(options['tableName']):
			table = group._f_getChild(options['tableName'])
			#print "table", options['tableName'], "already exists"
		else:
			table = self._h5file.createTable(group, options['tableName'], options['tableStruct'], options['tableTitle'], expectedrows=300000)
			#print "table", options['tableName'], "created"
		
		table.append(data)
		table.flush()

		#print "the data has been inserted into table", options['tableName']

		return table
		
	def addCategory(self, **groupinitparams):
		options = {
			groupName : 'Protein', 
			groupDesc : 'properties'
		}
		options.update(groupinitparams)
		
		assert self._h5file.root.__contains__(options['groupName']) == False, "group is already present, add another group with a different name"
		
		return self._h5file.createGroup(_h5file.root, option['groupName'], option['groupDesc'])


# This class appends new data into underlying tables
# does not yet implement /dev/shm

class Loader:
	def __init__(self, location, analysis_group='Protein'):
		# location is the execution path of this class
		self._location = location
		self._xvgfile = xvgfile.XVGFile()
		self._result = result.Result(location)
		self._analysis_group = analysis_group
		logging.basicConfig(filename='loader.log',level=logging.DEBUG)

	def load(self, filename, analysisName, TableStructure, fixed):
		# filename is the name of the analysis file
		# analysisName is the name of the analysis eg. 'rg'
		# TableStructure is the rowtype object to be created and inserted into the table
		# fixed number of fixed data columns in table (ie. info is not read in from 'filename'
		
		if analysisName.find(".") != -1:
			name,ext=analysisName.split(".")
		else:
			name = analysisName

		analysisRoot = os.path.join(self._location, name)
		xvgfilepath = os.path.join(analysisRoot, filename)
		
		logging.info("Loader.load: parsing %s", xvgfilepath)
		
		data = self._xvgfile.parse(TableStructure, fixed, xvgfilepath)
		table = self._result.addToTable(data, self._analysis_group, tableName=analysisName, tableStruct=TableStructure)

