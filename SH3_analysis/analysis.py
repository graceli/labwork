# -*- coding: utf-8 -*-                                                                                                                         
#!/usr/bin/env python
import tables 
import numpy
import glob
import sys
import os
import matplotlib.pylab as plt
import ConfigParser
import os

def parseConfig():
	#look for a configuration file in current directory
	#if not found then looking of ~/.analysisrc
	#if not found then quit with error
	#if no config in current dir then plot with default params
	config = ConfigParser.ConfigParser()
	default = '~/config.ini'
	current = os.path.join(os.getcwd(), 'config.ini')
	config.readfp(open(os.path.expanduser(default)))
	config.read(current)

	if not os.path.exists(current):
		print "you did not have a ini file in the current directory,",
	default, "will be used to set the entire analysis settings"

	return config

def plot_temp(plt, temp, na):
	avgIndex = 0
	stdIndex = 1
	y = na[:,avgIndex]
	err = na[:,stdIndex]
	plt.title(temp)
	plt.xlabel('sampling blocks')
	plt.errorbar(range(1,len(y)+1), y, yerr=err)


if len(sys.argv) < 3:
	print "usage: analysis.py path propertyname"
	sys.exit(0)

config = parseConfig()
path = sys.argv[1]
group, tablename, column = path.split('/')[1:]
print group, tablename, column, path
pname = sys.argv[2]

filelist = glob.glob("*.h5")

assert len(filelist) > 0

Tlist = config.get('system-specific', 'temps').split()
#numpy.genfromtxt('templist.txt')
print Tlist

count = 0
start = 1
end = 41

property = {}
for i in range(start, end):
	analysisfile = "analysis_part"+str(i)+".h5"

	if not os.path.exists(analysisfile):
		print analysisfile, "does not exist"
		continue

	print "processing", analysisfile
	file = tables.openFile(analysisfile)
	table = file.getNode('/'+group+'/'+tablename)
	data = []
	for Tstr in Tlist:
		T = int(Tstr)
		data = table.readWhere('temp==%(T)d' % vars())[column]
		if T not in property:
			property[T] = [[numpy.average(data), numpy.std(data)]]
		else:
			property[T].append([numpy.average(data), numpy.std(data)])

	count+=1

na270 = numpy.array(property[270])
na300 = numpy.array(property[300])
na450 = numpy.array(property[450])
na640 = numpy.array(property[640])
	
if config.get('outputting', 'txt'):
	plt.savetxt('na270.txt', na270, fmt='%0.2f %f')
	plt.savetxt('na300.txt', na300, fmt='%0.2f %f')
	plt.savetxt('na450.txt', na450, fmt='%0.2f %f')
	plt.savetxt('na640.txt', na640, fmt='%0.2f %f')


if config.get('outputting', 'plot'):
	#set plot settings
	if config.get('plotting', 'axes_tight'):
		plt.axis('tight')
	plt.grid(b=config.get('plotting', 'grid'), which='both')	
	plt.rcParams['savefig.dpi'] = config.get('plotting', 'resolution')
	plt.rcParams['font.size'] = config.get('plotting', 'font_size')

	#plot
	plt.subplot(221)
	plot_temp(plt, '270', na270)
	plt.savefig('testtest.png')
