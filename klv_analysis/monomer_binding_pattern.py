import csv
import sys
import glob
import math


# open X number of files for the monomer (there are 500)
# chiro_sys249_per_res_contacts.dat
fileList = glob.glob("*per_res_contacts.dat")
TOTAL_FILES = len(fileList)
#print fileList

#keep 6 arrays
#list of list

NRES = 7
segments = [ [0]*NRES, [0]*NRES , [0]*NRES, [0]*NRES, [0]*NRES]
average = [0]*NRES

seg = 0
filesProcessed = 0
for file in fileList:
	print "processing", file

	r = csv.reader(open(file), delimiter=' ')	
	nrows = 0
	for row in r:
		#print row
		time = float(row[0])
		for i in range(1,NRES+1): 
			value = int(row[i])
			segments[seg][i-1] += value
			average[i-1] += value 

		nrows += 1
	
	filesProcessed += 1
	print "processed", filesProcessed

	if filesProcessed == 100:
		seg += 1	
		print "at seg", seg
		filesProcessed = 0

#print seg
#sys.exit()
		
# compute the standard deviation for each value in the array 
sd = [0]*7
TOTAL_FRAMES = nrows
NSEGS = seg

print "total frames", TOTAL_FRAMES, "NSEGS", seg

for i in range(0, NRES):
	for s in range(0, seg):
		sd[i] += math.pow(float(segments[s][i])/(100*TOTAL_FRAMES) - float(average[i])/(TOTAL_FILES*TOTAL_FRAMES),2)
	sd[i] = math.sqrt(sd[i]/NSEGS)


# output the file 
# one column for the final values, second column is the stdev
print "# residue average sd"
print "# files Processed", filesProcessed
print "# total frames per file", TOTAL_FRAMES
print "$ total segments", seg
for i in range(0,NRES):
	print i, float(average[i])/(TOTAL_FILES*TOTAL_FRAMES), sd[i]

