import csv
import math
import sys

filename = sys.argv[1]

#constants
KB = 0.008314277
KB = KB*4.184
#r = csv.reader(open('AVGE_values.prn'), delimiter='\t')
f = open(filename)
tempList=[]
PElist=[]
for row in f:
	v = row.split()
	T=int(v[0])
	tempList.append(T)

	averagePE = float(v[1])/4.184
	PElist.append(averagePE)	

numTemps = len(tempList)
avalues = [0]*numTemps

sys.stderr.write("#temperature avalues (Kcal/mol)")
sys.stderr.write("printing the A values in the form of DR script file")

for i in range(1,numTemps):
	T = tempList[i]
	Tprev = tempList[i-1]
	beta_prev = 1/(KB*Tprev)
	beta = 1/(KB*T)
	
	avalues[i] = (beta_prev*avalues[i-1] + (beta-beta_prev)*0.5*(PElist[i-1]+PElist[i]))/beta
	#print T, avalues[i]

l = range(1,numTemps)
l.reverse()
for i in l:
	print "JOB", tempList[i], "5000 3", avalues[i]

print "JOB", tempList[0], "5000 3", 0.000

