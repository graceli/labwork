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
	
def run(f):
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
