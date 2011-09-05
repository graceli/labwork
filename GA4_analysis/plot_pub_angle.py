import pylab
import numpy
import analysis_config

def plot_bound_angles(fname, color, label):
	data = numpy.genfromtxt(fname+'.txt')
	x = data[:, 0]
	y = data[:, 1]
	result, bins = numpy.histogram(y, bins=30, normed=True)
	pylab.plot(bins[0:bins.shape[0]-1], result, color=color, label=label)
	
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
	for iso in ["scyllo", "chiro"]:
		fname = iso + '_data' + '_' + '0.35'
		color = analysis_config.LINE_COLOR[iso]
		plot_bound_angles(fname, color, iso)
		# plot_hist2d(fname)
	
	pylab.legend()
	pylab.grid(True)
	pylab.savefig('pub_angle_distribution.png', dpi=300)
	
	
if __name__ == '__main__':
	main()
