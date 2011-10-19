import numpy
import pylab
import glob

def plot_hist(name, files):
	fig = pylab.figure()
	ax = fig.add_subplot(111)
	data= []
	for f in files:
		data = numpy.genfromtxt(f, skiprows=12)
		# n, bins, patches = pylab.hist(data[:,1], bins=10, normed=1, align='left')
		hist, bins = numpy.histogram(data[:,1], bins=numpy.arange(0,1.5,0.05), normed=1)
		print hist
		print bins
		ax.plot(bins[0:len(hist)], hist)
		ax.set_xticks(bins[0::10])
		majorFormatter = pylab.FormatStrFormatter('%0.2f')
		minorLocator   = pylab.MultipleLocator(0.1)
		ax.xaxis.set_major_formatter(majorFormatter)
		ax.xaxis.set_minor_formatter(majorFormatter)
		ax.xaxis.set_minor_locator(minorLocator)
		
		ax.grid(True)
		# ax.set_xlabel('RMSD (nm)',size='8')
		
	fig.savefig(name+'.pdf')

def main():
	"""docstring for main"""

	files_15=glob.glob("*15*.gz")
	files_64=glob.glob("*64*.gz")

	pylab.rcParams['xtick.labelsize']='8'
	pylab.rcParams['ytick.labelsize']='8'
	pylab.rcParams['legend.fontsize']='8'	
	pylab.axis('tight')

	plot_hist('15', files_15)
	plot_hist('64', files_64)

if __name__ == '__main__':
	main()
	