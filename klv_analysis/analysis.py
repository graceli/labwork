#!/usr/bin/env python

import os
import sys
import pylab
import numpy
#how do I define my own import libraries?
#import "/work/grace/AnalysisScripts/pylibs/plotconfig" as plotconfig

def columnAverage(data,colnum):
	total_num_systems = len(data)
	stats = []
	ndatapoints = len(data[0][:,1])
	print "there are",ndatapoints, "datapoints" 
	for time in range(0, ndatapoints):
		values = []
		for i in range(0, total_num_systems):	
			if time < len(data[i][:,colnum]):
				values.append(data[i][time,colnum])

#		datapoint_time = data[i][time,0]	
		avg = numpy.average(values)
		std = numpy.std(values)
		stats.append([avg, std])	

	print "computed", len(stats), "number of points"

	return numpy.array(stats)	


#filebasename = sys.argv[1]
scyllodata = []
chirodata = []
#inositol_100mM_chiro_sys0_nosol.xtc_whole.xtc_p2p_vs_t.dat
for i in range(0, 10):
	scyllodatafile =  "inositol_100mM_scyllo_sys%(i)s_nosol.xtc_whole.xtc_p2p_vs_t.dat" % vars()
	if os.path.exists(scyllodatafile):
		scyllodata.append(numpy.genfromtxt(scyllodatafile))
	else:
		print "Error: did not find", scyllodatafile
	
	chirodatafile = "inositol_100mM_chiro_sys%(i)s_nosol.xtc_whole.xtc_p2p_vs_t.dat" % vars()
	if os.path.exists(chirodatafile):
		chirodata.append(numpy.genfromtxt(chirodatafile))
	else:
		print "Error: did not find", chirodatafile

inter_col = 1
intra_col = 2

scyllo_interhb = columnAverage(scyllodata, inter_col)
scyllo_intrahb = columnAverage(scyllodata, intra_col)
chiro_interhb = columnAverage(chirodata, inter_col)
chiro_intrahb = columnAverage(chirodata, intra_col)

pylab.subplot(221)
pylab.title("scyllo interhb vs t")
pylab.plot(scyllo_interhb[:,0])
numpy.savetxt('scyllo_interhb.txt', scyllo_interhb, fmt="%f %f")
pylab.subplot(222)
pylab.title("scyllo intrahb vs t")
pylab.plot(scyllo_intrahb[:,0])
numpy.savetxt('scyllo_intrahb.txt', scyllo_intrahb, fmt="%f %f")
pylab.subplot(223)
pylab.title("chiro interhb vs t")
pylab.plot(chiro_interhb[:,0])
numpy.savetxt('chiro_interhb.txt', chiro_interhb, fmt="%f %f")
pylab.subplot(224)
pylab.title("chiro intrahb vs t")
pylab.plot(chiro_intrahb[:,0])
numpy.savetxt('chiro_intrahb.txt', chiro_intrahb, fmt="%f %f")

pylab.savefig("so_interpeptide_hbond_ts.png")


