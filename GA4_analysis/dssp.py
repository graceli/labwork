import plot_and_save2hdf5 as myh5
import tables
# import numpy
import glob
import os

def process_dssp(filename, totalResidue, correction_factor, h5file='analysis_results.h5'):
	fp = open(filename)

	#initialize structure lists
	legend={}
	averageStruct = {}
	columnTotal = 0
	columnIndex = 0
	totalFramesProcessed=0
	raw_data = []
	for line in fp:
		if line[0] == "#":
			continue;
		elif line[0] == "@":
			columns = line.split()
			#print columns
			if columns[1][0] == "s" and columns[1] != "subtitle":
				#print columns
				structureType = columns[3][1:len(columns[3])-1]
				#print structureType
				legend[columnIndex+1] = structureType
				columnIndex+=1
				#print columnIndex

			columnTotal = columnIndex
			#initialize data array 
			for i in range(1, columnTotal+1):
				averageStruct[i]=0
		else:
			# should all be data now
			cols = line.split()
			raw_data.append(cols)
			for i in range(1,columnTotal+1):
				# correct for the 3 extra residues are counted in the GA4 system by dssp
				if legend[i] == "Coil":
					averageStruct[i] += (float(cols[i]) - correction_factor)/totalResidue
				else:
					averageStruct[i] += float(cols[i])/totalResidue
			totalFramesProcessed+=1

	# print "total number of columns is", columnTotal
	table = []
	table_descr = {}
	table.append(filename)
	table_descr['filename'] = tables.StringCol(256, pos=0)
	
	for i in range(1,columnTotal+1):
		table.append(averageStruct[i]/totalFramesProcessed)
		table_descr[legend[i]] = tables.Float32Col(pos=i)
	
	table.append(totalFramesProcessed)
	table_descr['num_frames'] = tables.Float32Col(pos=columnTotal+1)
	
	h5 = myh5.initialize(h5file)
	
	basename,ext = os.path.splitext(filename)
	myh5.save(h5, [tuple(table),], '/dssp/%(basename)s' % vars(), table_descr)
	
	raw_data_array = numpy.Array(raw_data)
	(nrows, ncols) = raw_data_array.shape
	myh5.save(h5, raw_data_array, '/dssp_data/%(basename)s' % vars(), myh5.create_description('col', ncols, format=tables.Int32Col(dflt=0)))
	

if __name__ == '__main__':
	filename='systematic_ap1f_scyllo_system4_scount.xvg'
	total_residues = 32
	correction_factor = 3
	files = glob.glob("*.xvg")
	for filename in files:
		print filename

		process_dssp(filename, total_residues, correction_factor)
	