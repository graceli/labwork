# -*- coding: utf-8 -*-
import numpy
import logging
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

		print "a h5 file", self._filename, "has been created in", self.location
		print "and a group with name",options['groupName'],"has been added"

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
			print "table", options['tableName'], "already exists"
		else:
			table = self._h5file.createTable(group, options['tableName'], options['tableStruct'], options['tableTitle'], expectedrows=300000)
			print "table", options['tableName'], "created"
	
		data_numpy = numpy.array(data)
		print data_numpy

		table.append(data_numpy)
		table.flush()

		print "the data has been inserted into table", options['tableName']

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
class Loader:
	def __init__(self, location, analysis_group='Protein'):
		# location is the execution path of this class
		self._location = location
		self._xvgfile = xvgfile.XVGFile()
		self._result = Result(location)
		self._analysis_group = analysis_group
		logging.basicConfig(filename='loader.log',level=logging.DEBUG)

	# def load(self, filename, analysisName, TableStructure, fixed):
	def load(self, analysis, xtc):
		# filename is the name of the analysis file
		# analysisName is the name of the analysis eg. 'rg'
		# TableStructure is the rowtype object to be created and inserted into the table
		# fixed number of fixed data columns in table (ie. info is not read in from 'filename'
		
		for filename, struct, analysis_name in zip(analysis.files(), analysis.types(), analysis.type_names()):
			xvgfilepath = os.path.join('/dev/shm', filename)
			logging.info("Loader.load: parsing and loading %s", xvgfilepath)
	
			print "parsing", filename	

			data = self._xvgfile.parse(struct, xvgfilepath, xtc)

			new_data_vector = xtc.params()	
			new_data_vector.append(data)
			# print new_data_vector
			table = self._result.addToTable(data, groupName=self._analysis_group, tableName=analysis_name, tableStruct=struct)







