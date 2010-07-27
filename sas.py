from tables import *
import csv
import os

class SAS(IsDescription):
	time = Int32Col() 
	replicanum = Int32Col()
	hydrophobic = Float32Col()
	hydrophilic = Float32Col()
	total = Float32Col()



def parseFilename(filename):
	stripped = os.path.basename(filename)
	parts = stripped.split("_")
	return [parts[0][4:], parts[1][7:]]

def addRecord(frame, replicanum, line):
	frame['time'] = int(line[0])	
	frame['replicanum'] = int(replicanum)
	frame['hydrophobic'] = float(line[1])
	frame['hydrophilic'] = float(line[2])
	frame['total'] = float(line[3])

def loadXVG2Table(table, partnum, replicanum, file):
	csvfile = open(file)
	r = csv.reader(csvfile, skipinitialspace=True, delimiter=' ')

	#parse out # and @ lines
	frame = table.row
	for line in r:
		if line[0][0] == "#" or line[0][0] == "@":
			print line
			continue
		else:
			print line
			addRecord(frame, replicanum, line)
			frame.append()

	table.flush()
	
					


filename = "sas.h5"

h5file = openFile(filename, mode = "w", title = "datafile")
group = h5file.createGroup("/", "Protein", "global properties")


# load the xvg data files
file="data/part1_replica101_sas.xvg"
partnum, replicanum = parseFilename(file)
table = h5file.createTable(group, 'SAS', SAS, "SAS data part %s by replica " % partnum)
loadXVG2Table(table, partnum, replicanum,file)

h5file.close()
#table = h5file.root.detector.readout
#test 

#example of ways to get stats
#numpy.histogram([row['hydrophobic'] for row in ro], bins=numpy.arange(23,36,0.1),normed=True)


loadSAS()
loadRadiusOfGyration()


