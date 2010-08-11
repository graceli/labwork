# -*- coding: utf-8 -*-                                                                                                                         
import tables 
import numpy
import glob
import sys
import os
import matplotlib.pylab as plt

if len(sys.argv) < 3:
	print "usage: analysis.py path propertyname"
	sys.exit(0)

path = sys.argv[1]
group, tablename, column = path.split('/')[1:]
print group, tablename, column, path
pname = sys.argv[2]

filelist = glob.glob("*.h5")

assert len(filelist) > 0

Tlist = numpy.genfromtxt('templist.txt')

count = 0

start = 1
end = 41
property = {}
#for analysisfile in filelist:
for i in range(start, end):
	analysisfile = "analysis_part"+str(i)+".h5"

	if not os.path.exists(analysisfile):
		print analysisfile, "does not exist"
		continue

	print "processing", analysisfile
	file = tables.openFile(analysisfile)
	table = file.getNode('/'+group+'/'+tablename)
	data = []
	for T in Tlist:
		data = table.readWhere('temp==%(T)d' % vars())[column]
		if T not in property:
			property[T] = [[numpy.average(data), numpy.std(data)]]
		else:
			property[T].append([numpy.average(data), numpy.std(data)])

	#count+=1

	#numpy.savetxt('%(pname)s_%(count)d.txt' % vars(), property, fmt='%d %f %f')

#print property

avgIndex = 0
stdIndex = 1

na270 = numpy.array(property[270])
#savetxt('na.txt', na, fmt='%0.1f %f')

plt.savetxt('na270.txt', na270, fmt='%0.2f %f')

na300 = numpy.array(property[300])
plt.savetxt('na300.txt', na300, fmt='%0.2f %f')

na450 = numpy.array(property[450])
plt.savetxt('na450.txt', na450, fmt='%0.2f %f')

na640 = numpy.array(property[640])
plt.savetxt('na640.txt', na640, fmt='%0.2f %f')

#plt.axis('tight')
plt.subplot(221)
y = na270[:,avgIndex]
err = na270[:,stdIndex]
title('270')
ylabel('sampling blocks')
plt.errorbar(range(1,len(y)+1), y, yerr=err)

plt.subplot(222)
y = na300[:,avgIndex]
err = na300[:,stdIndex]
title('300')
ylabel('sampling blocks')
plt.errorbar(range(1,len(y)+1), y, yerr=err)


#plt.subplot(223)
#y = na350[:,avgIndex]
#err = na350[:,stdIndex]
#plt.errorbar(range(1,len(y)+1), y, yerr=err)

plt.subplot(223)
y=na450[:,avgIndex]
err=na450[:,stdIndex]
plt.errorbar(range(1,len(y)+1), y, yerr=err)

plt.savefig(pname)

