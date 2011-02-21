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