import pylab
import numpy
import analysis_config

def plot_bound_angles(fname, color, line_type, label):
	data = numpy.genfromtxt(fname+'.txt')
	x = data[:, 0]
	y = data[:, 1]
	result, bins = numpy.histogram(y, bins=30, normed=True)
	
	lw = 1.0
	if line_type == '--':
		lw = 2.0
	
	pylab.plot(bins[0:bins.shape[0]-1], result, linewidth=lw, linestyle=line_type, color=color, label=label)
	
def plot_hist2d(fname):
	data = numpy.genfromtxt(fname + '.txt')
	x = data[:, 0]
	y = data[:, 1]
	H, xedges, yedges = numpy.histogram2d(x, y, bins=(30,30))
	print H
	extent = [yedges[0], yedges[-1], xedges[0], xedges[-1] ]
	pylab.imshow(H, extent=extent, aspect=1000)
	

def main():
	"""docstring for main"""
	analysis_config.plot_settings(setting="pub")
	# pylab.ylabel(r'P($\alpha$)')
	# pylab.xlabel(r'$\alpha(\deg)$')
	
	line_type = {'scyllo':'solid' , 'chiro': '--'}
	for iso in ["scyllo", "chiro"]:
		fname = iso + '_data' + '_' + '0.35'
		#color = analysis_config.LINE_COLOR[iso]
		lt = line_type[iso]
		plot_bound_angles(fname, 'black', lt, iso)
		# plot_hist2d(fname)
	
	pylab.legend()
	pylab.grid(True)
	pylab.savefig('pub_angle_distribution.pdf', dpi=1200)
	
	
if __name__ == '__main__':
	main()
