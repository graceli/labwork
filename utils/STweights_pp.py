import csv
import math

#constants
KB = 0.008314277
KB = KB*4.184
#r = csv.reader(open('AVGE_values.prn'), delimiter='\t')
f = open('AVGE_values.prn')
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
print "#temperature avalues (KJ/mol) avalues(Kcal/mol) averagePE(KJ/mol)"
print tempList[0], 0, 0, PElist[0]-PElist[0]
for i in range(1,numTemps):
	T = tempList[i]
	Tprev = tempList[i-1]
	beta_prev = 1/(KB*Tprev)
	beta = 1/(KB*T)
	
	avalues[i] = (beta_prev*avalues[i-1] + (beta-beta_prev)*0.5*(PElist[i-1]+PElist[i]))/beta
	print T, avalues[i], avalues[i]/4.184, PElist[i]-PElist[0]
