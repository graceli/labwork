import pylab
import os
import numpy
import tables
import plot_and_save2hdf5 as myh5

def plot_hist(ax, filename, label):
	if os.path.exists(filename):
		data = numpy.genfromtxt(filename)
		result, bins = numpy.histogram(data[:,2], bins=30, normed=True)
		ax.plot(bins[0:bins.shape[0]-1], result)
		ax.legend((label,), 'upper left', ncol=2)
	else:
		print filename, "does not exist"

def plot_dist(ax, filename, label):
	if os.path.exists(filename):
		data = numpy.genfromtxt(filename)
		#print data.shape
		ax.plot(data[:,0], data[:,1])
		ax.legend((label,), 'upper left', ncol=2)
	else:
		print filename, "does not exist"

def save_helper(h5file, data_array_list, paths_list):
	# for data_array, path in zip(data_array_list, paths_list):
		# table = myh5.getTable(h5file, path)
	for i in range(0, len(paths_list)):
		#print "data_array before saving"
		#print data_array_list[0]
		myh5.save(h5file, data_array_list[i], paths_list[i], format="float")

def read_helper(filename, big_array):
	results = numpy.array
	if os.path.exists(filename):
		results = numpy.genfromtxt(filename)
	else:
		print filename, "does not exist"
		return

	if big_array.shape[0] == 0:
		big_array = results
	else:
		big_array = numpy.append(big_array, results, axis=0)

	return big_array

def plot_from_file():
	""" 
		Plots angle analysis for GA4 beta oligomer systems
	"""

	pylab.rc('text', usetex=True)
	pylab.title(r'(GA)$_4$ $\beta$-sheet (all data)')
	pylab.xlabel(r'angle ($\deg$)')
	pylab.ylabel('Population (normalized by total)')

	isomer = ["scyllo", "chiro"]
	analysis = ["angle", "dist"]
	
	angle_fig = pylab.figure(); angle_ax = angle_fig.add_subplot(111)
	dist_fig = pylab.figure(); dist_ax  = dist_fig.add_subplot(111)
	
	for iso in isomer:
		for i in range(0, 3):
			angle_name = 'GA4_perfect_' + iso + str(i) + '_nosol.xtc_angle.xvg'
			dist_name = 'GA4_perfect_' + iso + str(i) + '_nosol.xtc_dist.xvg'
			
			plot_hist(angle_ax, angle_name, iso + str(i))
			plot_dist(dist_ax, dist_name, iso + str(i))
	
	angle_fig.savefig('inositol_angle_distr_for_GA4_perfect.png')
	savefig('inositol_dist_for_GA4_perfect.png')

def read_and_save_data():
	h5file = tables.openFile('GA4_angles_analysis.h5', mode='a')
	for iso in ["scyllo", "chiro"]:
		for sys in range(0,3):
			angle_all = numpy.array([])
			mindist_all = numpy.array([])
			for i in range(0,4):			
				# put all inositol analysis in a single table per system
				# have separate tables for angles and distances_from_sheet
				angle_name = 'GA4_perfect_' + iso + str(sys) + '_nosol.xtc_angle' + '_inos' + str(i) + '.xvg'
				mindist_name = 'GA4_perfect_' + iso + str(sys) + '_nosol.xtc_mindist' + '_inos' + str(i+14) + '.xvg'
				
				angle_all = read_helper(angle_name, angle_all)				
				mindist_all = read_helper(mindist_name, mindist_all)
				
				if angle_all == None or mindist_all == None:
					break

			if angle_all is not None and mindist_all is not None and angle_all.shape[0] > 0 and mindist_all.shape[0] > 0:
				paths = [ '/angle/%(iso)s%(sys)d' % vars(), '/mindist/%(iso)s%(sys)d' % vars() ]
				save_helper(h5file, [angle_all, mindist_all], paths)
	return h5file
	
def get_angle_matrix(h5file, sys):
	"""docstring for get_angle"""
	path = os.path.join('/angle', sys)
	print "getting ", path
	matrix = myh5.getTableAsMatrix(h5file, path, dtype=numpy.float64)
	return matrix
	
def get_mindist_matrix(h5file, sys):
	"""docstring for get_mindist_matrix"""
	path = os.path.join('/mindist', sys)
	print "getting ", path
	
	return myh5.getTableAsMatrix(h5file, path, dtype=numpy.float64)

def plot_bound_angles(fname):
	data = numpy.genfromtxt(fname+'.txt')
	x = data[:, 0]
	y = data[:, 1]
	result, bins = numpy.histogram(y, bins=30, normed=True)
	fig = pylab.figure()
	ax = fig.add_subplot(111)

	ax.plot(bins[0:bins.shape[0]-1], result)
	ax.legend((r'P($\alpha$)',), 'upper left', ncol=2)
	fig.savefig(fname + '_angle_distribution.png')

def plot_hist2d(fname):
	data = numpy.genfromtxt(fname + '.txt')
	x = data[:, 0]
	y = data[:, 1]
	H, xedges, yedges = numpy.histogram2d(x, y, bins=(30,30))
	print H
	extent = [yedges[0], yedges[-1], xedges[0], xedges[-1] ]
	pylab.imshow(H, extent=extent, aspect=1000)
	pylab.savefig(fname + "_image.png")


def analysis(h5file):
	print "performing analysis"
	cutoff = 0.35
	for iso in ["scyllo", "chiro"]:
		fname = iso + '_data' + '_' + str(cutoff)
		f = open(fname + '.txt', 'w')
		for s in range(0,3):
			angle_mat = get_angle_matrix(h5file, iso+str(s))
			mindist_mat = get_mindist_matrix(h5file, iso+str(s))
			# smarter: use join where mindist_col > 0.5
			# get only angles where mindist to fibril is less than 0.6
			if angle_mat is not None or mindist_mat is not None:
				nrows = min(angle_mat.shape[0], mindist_mat.shape[0])
				print angle_mat
				print mindist_mat	
				for i in range(0, nrows):
					angle = angle_mat[i][2]
					mindist = mindist_mat[i][1]
					if mindist < cutoff:
						#correct for angles
						if angle > 90:
							print >> f, mindist, angle - 90
						else:
							print >> f, mindist, angle

		f.flush()
		plot_bound_angles(fname)
		plot_hist2d(fname)

def main():
	"""docstring for main"""
	# h5file = read_and_save_data()
	h5file = tables.openFile('GA4_angles_analysis.h5')
	analysis(h5file)
	# for iso in ["scyllo", "chiro"]:
	# 	fname = iso + '_data' + '_' + '0.35'
	# 	plot_hist2d(fname)
	
if __name__ == '__main__':
	"""	Analysis:
		- Obj of analysis: Calculate planar angle between plane of inositol and plane of sheet.
		- Purpose: To understand how binding geometry differs between stereisomers to the sheet and whether this depends on the type of contact made.
	"""
	
	main()