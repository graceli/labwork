#!/usr/bin/env python

import sys
import math

#histogram

if len(sys.argv) < 4:
	print "usage <file> <binwidth> <norm> <time|notime>"
	exit(0)


fp = open(sys.argv[1])
binwidth = float(sys.argv[2])
norm = sys.argv[3]
maxBin=0
minBin=10000
histogram = {}
linenum=0

for line in fp:

        if line.find("#") == 0:
                continue

	if line.find("@") == 0:
		continue

	columns = line.split()
	numCols = len(columns)

	if numCols == 1:
		datapoint = float(columns[0])
	else:
		datapoint = float(columns[1])

	#map data to bin number
	binnum = int(math.floor(datapoint/binwidth))
	
	#determine max and min bin
	if binnum > maxBin:
		maxBin = binnum
	
	if binnum < minBin:
		minBin = binnum
	
	if binnum in histogram:
		histogram[binnum]+=1
	else:
		histogram[binnum]=1

	linenum+=1

if norm == "-norm":
	doNormalize = linenum
else:
	doNormalize = 1

print "# Detected", numCols, "but only calculating on column 1"
print "# binwidth", binwidth
print "# norm ", norm
#print "# unit ", unit
for i in range(minBin, maxBin):
        binnum = i+1
        datacenter = (binnum*binwidth+(binnum-1)*binwidth)/2
        if i in histogram:
                print datacenter, float(histogram[i])/doNormalize        
        else:
                print datacenter, 0.0
