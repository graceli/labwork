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

analysisParams = {'txtext':'csv', 'figext':'png'}

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
	
def setConfig(config):
	global analysisParams
	
	if config.get('outputting', 'txt'):
		txtformat = config.get('formatting', 'text')
		analysisParams['txtext'] = txtformat
		
	if config.get('outputting', 'plot'):
		#set plot settings
		if config.get('plotting', 'axes_tight'):
			plt.axis('tight')

		plt.grid(b=config.get('plotting', 'grid'), which='both')	
		plt.rcParams['savefig.dpi'] = config.get('plotting', 'resolution')
		plt.rcParams['font.size'] = config.get('plotting', 'font_size')
		figformat = config.get('formatting', 'figure')
		analysisParams['figext'] = figformat

def plot_by_parts(na270, na300, na450, na640):
	if config.get('outputting', 'txt'):
		plt.savetxt('na270.' + analysisParams['txtext'], na270, fmt='%0.6f %f')
		plt.savetxt('na300.' + analysisParams['txtext'], na300, fmt='%0.6f %f')
		plt.savetxt('na450.' + analysisParams['txtext'], na450, fmt='%0.6f %f')
		plt.savetxt('na640.' + analysisParams['txtext'], na640, fmt='%0.6f %f')

	if config.get('outputting', 'plot'):
		#plot hardcoded settings
		plt.subplots_adjust(top=0.95,hspace=0.4)

		plt.subplot(221)
		plot_vector(plt, '270', na270)
		plt.subplot(222)
		plot_vector(plt, '300', na300)
		plt.subplot(223)
		plot_vector(plt, '450', na450)
		plt.subplot(224)
		plot_vector(plt, '640', na640, 'sampling blocks')
		plt.savefig(pname + '.' + analysisParams['figext'])

def plot_hist(data):
	if config.get('outputting', 'plot'):
		pdf, bins, patches = plt.hist(data[300], arange(data[300].min(),
data[300].max(), 0.05), normed=True, facecolor='green', alpha=0.75)
	#n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)
	
	if config.get('outputting', 'txt'):
		s= numpy.vstack((pdf,bins[0:len(bins)-1])).transpose()
		plt.savetxt(pname+'_distribution' + analysisParams['txtext'], s,
fmt='%f %f')

def plot_vector(plt, title, na, label=''):
	#plt is the plot object
	#title is the plot title
	#na is the numpy array to be plotted
	
	avgIndex = 0
	stdIndex = 1
	y = na[:,avgIndex]
	err = na[:,stdIndex]
	plt.title(title)
	plt.errorbar(range(1,len(y)+1), y, yerr=err)
	plt.xlabel(label)
	
def plot_by_temp(by_temp):
	plt.savetxt('by_temp.csv', by_temp, fmt='%f %f')
	plot_vector(plt, 'Rg(Temperature)', by_temp, 'Temperature')
	plt.savefig('by_temp.png')
	
if __name__ == "__main__":
	if len(sys.argv) < 3:
		print "usage: analysis.py path propertyname"
		sys.exit(0)

	config = parseConfig()
	setConfig(config)

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

	by_parts = {}
	data = {}
	for i in range(start, end):
		analysisfile = "analysis_part"+str(i)+".h5"

		if not os.path.exists(analysisfile):
			print analysisfile, "does not exist"
			continue

		print "processing", analysisfile
		file = tables.openFile(analysisfile)
		atable = file.getNode('/'+group+'/'+tablename)
		
		for Tstr in Tlist:
			T = int(Tstr)
			datachunk = atable.readWhere('temp==%(T)d' % vars())[column]
			
			# add to or extend the current list with the data read from a table
			if T not in data:
				data[T] = datachunk
			else:
				data[T].extend(datachunk)
			
			# do per analysis file aggregation
			if T not in by_parts:
				by_parts[T] = [[numpy.average(datachunk), numpy.std(datachunk)]]
			else:
				by_parts[T].append([numpy.average(datachunk), \
				numpy.std(datachunk)])

		
		#close the h5 file
		file.close()

	by_temp = {}
	for Tstr in Tlist:
		T = int(Tstr)
		if T not in by_temp:
			by_temp[T] = [[numpy.average(data[T]), numpy.std(data[T])]]
		else:
			by_temp[T].append([numpy.average(data[T]), numpy.std(data[T])])

	na270 = numpy.array(by_parts[270])
	na300 = numpy.array(by_parts[300])
	na450 = numpy.array(by_parts[450])
	na640 = numpy.array(by_parts[640])

	plot_by_parts(na270,na300,na450,na640)
	plot_hist(data)
	plot_by_temp(by_temp)




