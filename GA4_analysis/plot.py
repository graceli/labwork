#!/usr/bin/env python
import pylab
import numpy
import matplotlib

def plot_settings(setting="default"):
	if setting == "pub":
		pylab.rcParams['figure.figsize'] = [6,3]
		pylab.rcParams['font.size'] = 9
		pylab.rcParams['legend.fontsize'] = 'small'
	else:
		pylab.rcParams['figure.figsize'] = [8,6]
		pylab.rcParams['font.size'] = 12


def plot_with_error(x, y, yerror, label, axes=None):
	""" plots a figure onto a figure object and 
	returns the modified figure"""
	
	axes.fill_between(x, y-yerror, y+yerror, alpha=0.1)
	axes.plot(x,y,label=label)
	# fig.savefig(name)
	
	return axes

def plot_eed(file_name="test.png"):
	# I don't have Latex installed
	# pylab.rcParams['text.usetex']=True
	# pylab.rcParams['text.latex.unicode']=True
	
	scyllo = numpy.genfromtxt('scyllo_eed_data2.dat', dtype=float)	
	chiro = numpy.genfromtxt('chiro_eed_data2.dat', dtype=float)	
	plot_settings(setting="pub")
	
	fig = pylab.figure()
	fig.subplots_adjust(top=0.9, bottom=0.19, left=0.18, right=0.90, wspace=0.001)

	ax1 = fig.add_subplot(1,2,1)
	# figure sharing doesn't work -- don't know why ...
	ax2 = fig.add_subplot(1,2,2)	
	ax1.set_xlabel(r'End to end distance ($\AA$)')
	ax1.set_ylabel("Probability")
	yticklabels = ax2.get_yticklabels()+ax2.get_xticklabels()
	pylab.setp(yticklabels, visible=False)

	nax = 0
	for eed_data in [scyllo, chiro]:
		nrows, ncols = eed_data.shape
		ax = fig.get_axes()[nax]
		ax.grid('on')
		ax.set_xlim(0, 25)
		ax.set_ylim(0, 0.1)
		for i in range(1, ncols, 2):
			print 0, i, i+1
			x = eed_data[:,0]
			y = eed_data[:,i]
			yerror = eed_data[:,i+1]
			state = int(i/2)
			ax = plot_with_error(x, y, yerror, '%(state)d hbonds' % vars(), axes=ax)
		nax += 1
	ax2.legend()
	fig.savefig("all.png")


def main():
	plot_eed()


if __name__ == '__main__':
	main()
