import plot_and_save2hdf5 as myh5
import tables
# import numpy
import glob
import os

def process_dssp(filename, totalResidue, h5file='analysis_results.h5'):
	fp = open(filename)

	#initialize structure lists
	legend={}
	averageStruct = {}
	columnTotal = 0
	columnIndex = 0
	totalFramesProcessed=0

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
			for i in range(1,columnTotal+1):
				# correct for the 3 extra residues are counted in the GA4 system by dssp
				if legend[i] == "Coil":
					averageStruct[i] += (float(cols[i]) - 3)/totalResidue
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
		table_descr[legend[i]] = tables.Float64Col(pos=i)
	
	table.append(totalFramesProcessed)
	table_descr['num_frames'] = tables.Float64Col(pos=columnTotal+1)
	
	h5 = myh5.initialize(h5file)
	
	# def save(h5file, data, table_path, table_struct=numpy.dtype(numpy.int64)):
	basename,ext = os.path.splitext(filename)
	myh5.save(h5, [tuple(table),], '/dssp/%(basename)s' % vars(), table_descr)
	

if __name__ == '__main__':
	filename='systematic_ap1f_scyllo_system4_scount.xvg'
	total_residues = 32
	files = glob.glob("*.xvg")
	for filename in files:
		print filename

		process_dssp(filename, total_residues)
	