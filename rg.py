import tables 
import csv 
import os

class RGTable(tables.IsDescription):
	time = tables.Int32Col()
	replicanum = tables.Int32Col()
	partnum = tables.Int32Col()
	rgval = tables.Float32Col()


def initTable(file, tablename, tabletitle, groupname, grouptitle, TABLETYPE):
	# assumes that this appends ??
	# can pass classes?
	table = file.createTable(file.root._f_getChild(groupname), tablename, TABLETYPE)
	
	return table


def loadIntoTable(table, partnum, replicanum, fileToLoad):
	csvfile = open(fileToLoad)
	r = csv.reader(csvfile, skipinitialspace=True, delimiter=' ')
	record = table.row
	for line in r:
		if line[0][0] == "#" or line[0][0] == "@":
			print line
			continue
		else:
			#print line
			setRgRecord(record, partnum, replicanum, line)
			record.append()

	table.flush()
	

#utility functions
def setRgRecord(record, partnum, replicanum, line):
	#TODO: should check for size of line and number entered into table
	
	record['time'] = float(line[0])	
	record['partnum'] = int(partnum)
	record['replicanum'] = int(replicanum)
	record['rgval'] = float(line[1])
		
def parseFilename(filename):
	stripped = os.path.basename(filename)
	name,ext = os.path.splitext(stripped)
	parts = name.split("_")
	return [parts[0][4:], parts[2]]


filename = "rg.h5"
tabletitle = "rg table"
groupname = "Protein"
grouptitle= ''
h5file = tables.openFile(filename, mode="w", title=tabletitle)
group = h5file.createGroup("/", groupname, grouptitle)	

dirList = os.listdir('rg')
tableList = []
for d in dirList:
	path = os.path.join('rg', d)
	filesList = os.listdir(path)
	table = initTable(h5file, "Rg_"+d, "rg data", "Protein", "global properties", RGTable)
	for filename in filesList:
		partnum,replicanum = parseFilename(filename)
		loadIntoTable(table, partnum, replicanum, os.path.join(path,filename))
		tableList.append(table)
		
h5file.close()

#open up and and load in all the rg datafile
#for d in dirList:
#	print d
#	filesList = os.listdir(d)
#	for filename in filesList:
#		print filename
		#loadIntoTable(table,filename)












