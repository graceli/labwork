import re
import sys

if len(sys.argv) < 3:
	print "use:", sys.argv[0], "<file input>", "<Nresidues>"
	exit(0) 
filename = sys.argv[1]
fp = open(filename)

#configurable variables
totalResidue = int(sys.argv[2])

#initialize structure lists
legend={}
averageStruct = {}
columnTotal = 0
columnIndex = 0
totalFramesProcessed=0

for line in fp:
	if line[0] == "#":
		continue;
	elif line[0] == "@":
		columns = line.split()
		#print columns
		if columns[1][0] == "s" and columns[1] != "subtitle":
			#print columns
			structureType = columns[3][1:len(columns[3])-1]
			#print structureType
			legend[columnIndex+1] = structureType
			columnIndex+=1
			#print columnIndex

		columnTotal = columnIndex
		#initialize data array 
		for i in range(1, columnTotal+1):
			averageStruct[i]=0
	else:
	# should all be data now
		cols = line.split()
		#print line
		#print cols[0]
		for i in range(1,columnTotal+1):
			averageStruct[i] += float(cols[i])/totalResidue
		totalFramesProcessed+=1

# print "total number of columns is", columnTotal
print filename
for i in range(1,columnTotal+1):
	print legend[i], averageStruct[i]/totalFramesProcessed, totalFramesProcessed
					
