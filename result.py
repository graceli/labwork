import tables
import rowtypes
import os

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
		
		#member variable naming conventions?
		self._filename = options['filename'] 
		self.targetpath = targetpath
		self.location = os.path.join(targetpath, self._filename)
		self.templocation = os.path.join('/dev/shm', self._filename)
		self._h5file = tables.openFile(self.templocation, mode = "a", title=options['tabletitle'])
		
		if not self._h5file.root.__contains__(options['groupName']):
			self._h5file.createGroup(self._h5file.root, options['groupName'], options['groupDesc'])

		print "a h5 file", self._filename, "has been created in", self.templocation
		print "and a group with name",options['groupName'],"has been added"

	def __del__(self):
		self._h5file.close()	
		os.system("cp %(self.templocation)s %(self.location)s" % vars())
		print "the file", self._filename,"has been copied to", self.location
		
	def addToTable(self, data, **tableinitparams):
		#dictionary embedded within a dictionary?
		options = {
				'groupName' : 'Protein',
				'tableName' : 'defaultname',
				'tableTitle' : 'defaulttitle',
				'tableStruct' : rowtypes.DefaultTable
		}
		options.update(tableinitparams)
		print "group to be added to is", options['groupName']
	
		assert self._h5file.root.__contains__(options['groupName']) == True, "group is not found add the group first"
			
		group = self._h5file.root._f_getChild(options['groupName'])

		if group.__contains__(options['tableName']):
			table = group._f_getChild(options['tableName'])
			print "table", options['tableName'], "already exists"
		else:
			table = self._h5file.createTable(group, options['tableName'], options['tableStruct'], options['tableTitle'])
			print "table", options['tableName'], "created"
		
		table.append(data)
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

if __name__ == "__main__":
	aresult = Result()		
	dummyrows = [(0,1,1000,300),(1,1,1001,320)]
	aresult.addToTable([(0,1,1000,300)])
	aresult.addToTable(dummyrows)
