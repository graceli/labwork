#!/usr/bin/env python
import os
import sys
import re
import csv
import glob
import pylab
import numpy
import tables
import plot_and_save2hdf5 as myh5
import utils

# iterate over all table objects
# for node in h5file.walkNodes("/"):
# 	print node #table object
 	# node._v_parent._v_name # name of the group

# columnAverage is now in utils.py

# inositol_so_cluster.py
def average_cluster_columns(datafile, outputname):
    """ this function assumes that every odd column starting from 1 is time or is not used"""
    data = numpy.genfromtxt(datafile)
    sdata = data[0:140001,1::2]
    print sdata
    averaged_columns = numpy.average(sdata, axis=1)
    print averaged_columns
    numpy.savetxt(outputname, numpy.transpose([data[0:140001,0],averaged_columns]), fmt='%0.2f')

# average_cluster_columns('scyllo_45to4_nclust_140ns.dat', 'scyllo_45to4_nclust_140ns.txt')
# average_cluster_columns('chiro_45to4_nclust_140ns.dat', 'chiro_45to4_nclust_140ns.txt')

def intersection(h5file, ratio):
	"""docstring for intersection"""

	#nasty fix for different table names
	polarName = {'15to4' : 'whole_nosol_0-200ns_inos_total.dat', '45to4' : 'whole_nosol_0-200_inos_total.dat'}
	nonpolarName = {'15to4' : 'whole_nosol_0-200ns_per_inositol_contacts.dat', '45to4' : 'whole_nosol_0-200_per_inositol_contacts.dat'}
	
	isomerList = ["scyllo", "chiro"]
	polarPath = "/polar"
	nonpolarPath = "/nonpolar_residue"
	dataList = [['isomer', 'system#', 'polar_only', 'polar_and_nonpolar', 'nonpolar_only', 'total']]
	
	resultsWriter = csv.writer(open(ratio + '_intersection.txt', 'wb'), delimiter=' ')
	
	for iso in isomerList:
		data = []
		for sys in range(1,21):
			polarFile = os.path.join(polarPath, "%(iso)s_sys%(sys)s_%(ratio)s_" % vars() + polarName[ratio])
			polarMatrix = myh5.getTableAsMatrix(h5file, polarFile)

			nonpolarFile = os.path.join(nonpolarPath, "%(iso)s_sys%(sys)s_%(ratio)s_" % vars() + nonpolarName[ratio])
			nonpolarMatrix = myh5.getTableAsMatrix(h5file, nonpolarFile)
			
			if polarMatrix != None and nonpolarMatrix != None:
				rows, cols = nonpolarMatrix.shape
				polar_and_nonpolar = 0.0
				polar_only = 0.0
				nonpolar_only = 0.0
				for i in range(20001, rows):
					for j in range(1, cols):
						if polarMatrix[i][j] and nonpolarMatrix[i][j]:
							polar_and_nonpolar += 1
						elif polarMatrix[i][j]:
							polar_only += 1
						elif nonpolarMatrix[i][j]:
							nonpolar_only += 1

				total = polar_only + polar_and_nonpolar + nonpolar_only	
				dataList.append([iso, sys, polar_only / total, polar_and_nonpolar / total, nonpolar_only / total, total])
				data.append([polar_only / total, polar_and_nonpolar / total, nonpolar_only / total, total])
		
		print data
		print numpy.array(data)
		average = numpy.average(numpy.array(data), axis=0)
		std = numpy.std(numpy.array(data), axis=0)
		print average.tolist()
		print std.tolist()
		# 
		addToListAvg = [iso+' avg', 'all']
		addToListAvg.extend(average.tolist())
		addToListStd = [iso+' std', 'all']
		addToListStd.extend(std.tolist())
		
		dataList.append(addToListAvg)
		dataList.append(addToListStd)

	# numpy.savetxt("15to4_intersection.gz", dataList, fmt='%s %d %0.3f %0.3f %0.3f %d')
	resultsWriter.writerows(dataList)

def nonpolar_residue(h5file, tag):
	scyllo_pattern = re.compile(r'scyllo')
	chiro_pattern = re.compile(r'chiro')
	atype_pattern = re.compile(r'residue_contact')
	data_list = {'scyllo':[], 'chiro':[]}
	
	#fix this number for now
	N_datapoints = 98000
	for table in h5file.listNodes("/nonpolar_residue", 'Table'):
		if atype_pattern.search(table.name):
			table_path = os.path.join("/nonpolar_residue", table.name)
			data = myh5.getTableAsMatrix(h5file, table_path)

			sum_over_time = numpy.average(data[20000:N_datapoints, 1:], axis = 0)
			# print sum_over_time.size
			
			sum_over_time.shape = (sum_over_time.size / 4, 4)
			sum_over_peptides = numpy.sum(sum_over_time, axis = 1)
				
			if scyllo_pattern.search(table.name):
				data_list['scyllo'].append(sum_over_peptides)
			elif chiro_pattern.search(table.name):
				data_list['chiro'].append(sum_over_peptides)
			else:
				print "No pattern matches", table.name
	
	# save results to flat files
	for isomer in data_list.keys():
		nparray = numpy.array(data_list[isomer])
				
		# dump the list of results for each system
		numpy.savetxt('%(tag)s_nonpolar_residue_inositol_contact_%(isomer)s_counts.txt' % vars(), nparray, fmt='%0.8f')
		print "saved", isomer, "analysis with shape", nparray.shape

		# average over all the systems; each system is a row in nparray
		average = numpy.average(nparray, axis=0)
		std = numpy.std(nparray, axis=0)

		#save the normalized average and std
		numpy.savetxt('%(tag)s_nonpolar_residue_inositol_contact_%(isomer)s_avg_std.txt' % vars(), [average, std], fmt='%0.8f')
		

def dssp_helper(f, iso, system):
	# structList = ['Coil', 'B-Sheet', 'B-Bridge', 'Bend', 'Turn', 'A-Helix', '3-Helix']
	heading = ['#filename', 'Coil', 'B-Sheet', 'B-Bridge', 'Bend', 'Turn', 'A-Helix', '3-Helix']
	ss_distribution = {'#filename': '', 'Coil' : 0.0, 'B-Sheet' : 0.0, 'B-Bridge' : 0.0, 'Bend' : 0.0, 'Turn' : 0.0, 'A-Helix' : 0.0, '3-Helix' : 0.0}

	# use regex to match the tables for 'system' and 'iso' -- another better way? 
	pattern = re.compile(r"%(iso)s.*%(system)s" % vars())

	#write data to ascii csv to import into Excel
	dictFilename = '%(system)s_%(iso)s_dict.csv' % vars()
	fhandle = open(dictFilename, 'wb')
 	writer = csv.DictWriter(fhandle, heading, delimiter=" ")
	writer.writerow(dict(zip(heading, heading)))

	num_matches = 0
	for t in f.listNodes(where='/dssp_results'):
		filename = t.col('filename')[0]
		match = pattern.search(filename)
		if match:
			ss_distribution['#filename'] = filename
			print "processing", filename
			# compute the sum of ss distribution over
			# all independent runs
			for struct in ss_distribution.keys():
				if struct != '#filename':
					try:
						ss = t.col(struct)[0]
					except KeyError:
						print filename, struct,"does not exist in table"
					else:
						# print filename, struct, ss
						ss_distribution[struct] = float(ss)
				else:
					continue
			writer.writerow(ss_distribution)
			fhandle.flush()
	fhandle.close()
	#regenerate csv file and compute the average and std (over all independent systems)
	#TODO: figure out how to do this without read off of disk
	#get rid of the first column because its all strings (read in as NaN by Numpy)
	
def dssp(f):
	# "Coil": Float64Col(shape=(), dflt=0.0, pos=2),
	#   "B-Sheet": Float64Col(shape=(), dflt=0.0, pos=3),
	#   "B-Bridge": Float64Col(shape=(), dflt=0.0, pos=4),
	#   "Bend": Float64Col(shape=(), dflt=0.0, pos=5),
	#   "Turn": Float64Col(shape=(), dflt=0.0, pos=6),
	#   "A-Helix": Float64Col(shape=(), dflt=0.0, pos=7),
	#   "3-Helix": Float64Col(shape=(), dflt=0.0, pos=8),

	""" 
		reads the h5 file containing the secondary structure analysis timeseries
		and outputs a distribution for each datafile
		write the computed distribution to a csv file with headings
	"""

	for iso in ['scyllo', 'chiro']:
		for system in ['15to4', '45to4']:
			dssp_helper(f, iso, system)
	
	dssp_helper(f, 'water', "")

	heading = ['#filename', 'Coil', 'B-Sheet', 'B-Bridge', 'Bend', 'Turn', 'A-Helix', '3-Helix']
	ss_distribution = {'#filename': '', 'Coil' : 0.0, 'B-Sheet' : 0.0, 'B-Bridge' : 0.0, 'Bend' : 0.0, 'Turn' : 0.0, 'A-Helix' : 0.0, '3-Helix' : 0.0}
	
	#write to each csv file the average and standard deviation
	data = []
	for dictFilename in glob.glob("*.csv"):
		try:
			data = numpy.genfromtxt(dictFilename, comments='#')[:,1:]
		except IOError:
			print data
			print dictFilename, "size", os.path.getsize(dictFilename), "had a problem"


		# print data
		average = numpy.average(data, axis=0)
		print average
		stdev = numpy.std(data, axis=0)
		writer = csv.DictWriter(open(dictFilename, 'a'), heading, delimiter=" ")
		#append the average and stdev back into the csv file
		aggregate_data = {'average':average, 'stdev':stdev}
		for newhead in ['average', 'stdev']:
			index = 0
			for attr in heading:
				if attr == "#filename":
					ss_distribution[attr] = newhead
				else:
					print index
					ss_distribution[attr] = aggregate_data[newhead][index]
					index += 1
			writer.writerow(ss_distribution)


def main():
	""" """
	# print sys.argv[1]
	filename = sys.argv[1]
	h5file = tables.openFile(filename)
	ratio,ext = os.path.splitext(filename)
	# nonpolar_residue(h5file, ratio)
	# intersection(h5file, ratio)
	dssp(h5file)

if __name__ == '__main__':
	main()
	
# def main():
# 	"""docstring for main"""
# 	#filebasename = sys.argv[1]
# 	scyllodata = []
# 	chirodata = []
# 	#inositol_100mM_chiro_sys0_nosol.xtc_whole.xtc_p2p_vs_t.dat
# 	for i in range(0, 10):
# 		scyllodatafile =  "inositol_100mM_scyllo_sys%(i)s_nosol.xtc_whole.xtc_p2p_vs_t.dat" % vars()
# 		if os.path.exists(scyllodatafile):
# 			scyllodata.append(numpy.genfromtxt(scyllodatafile))
# 		else:
# 			print "Error: did not find", scyllodatafile
# 	
# 		chirodatafile = "inositol_100mM_chiro_sys%(i)s_nosol.xtc_whole.xtc_p2p_vs_t.dat" % vars()
# 		if os.path.exists(chirodatafile):
# 			chirodata.append(numpy.genfromtxt(chirodatafile))
# 		else:
# 			print "Error: did not find", chirodatafile
# 
# 	inter_col = 1
# 	intra_col = 2
# 
# 	scyllo_interhb = columnAverage(scyllodata, inter_col)
# 	scyllo_intrahb = columnAverage(scyllodata, intra_col)
# 	chiro_interhb = columnAverage(chirodata, inter_col)
# 	chiro_intrahb = columnAverage(chirodata, intra_col)
# 
# 	pylab.subplot(221)
# 	pylab.title("scyllo interhb vs t")
# 	pylab.plot(scyllo_interhb[:,0])
# 	numpy.savetxt('scyllo_interhb.txt', scyllo_interhb, fmt="%f %f")
# 	pylab.subplot(222)
# 	pylab.title("scyllo intrahb vs t")
# 	pylab.plot(scyllo_intrahb[:,0])
# 	numpy.savetxt('scyllo_intrahb.txt', scyllo_intrahb, fmt="%f %f")
# 	pylab.subplot(223)
# 	pylab.title("chiro interhb vs t")
# 	pylab.plot(chiro_interhb[:,0])
# 	numpy.savetxt('chiro_interhb.txt', chiro_interhb, fmt="%f %f")
# 	pylab.subplot(224)
# 	pylab.title("chiro intrahb vs t")
# 	pylab.plot(chiro_intrahb[:,0])
# 	numpy.savetxt('chiro_intrahb.txt', chiro_intrahb, fmt="%f %f")
# 
# 	pylab.savefig("so_interpeptide_hbond_ts.png")
# 
# 
# if __name__ == '__main__':
# 	main()

