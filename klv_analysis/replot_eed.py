import pylab
import numpy
import math
import config

def plot_with_error(x, y, yerror, label, color, facecolor, axes=None):
	""" plots a figure onto a figure object and 
	returns the modified figure"""
	
	axes.fill_between(x, y-yerror, y+yerror, alpha=0.5, facecolor=facecolor, edgecolor=facecolor)
	axes.plot(x, y, label=label, color=color)
	# fig.savefig(name)
	
	return axes

def plot_eed(file_name="test.png"):
	# plots the end to end distance distributions
	# used for the KLVFFAE peptide monomer system
	# flat file data in ~/Desktop/KLV_paper/data/analysis/monomer
	pylab.clf()
	params = {'text.usetex': False,  # turn this off to use Arial font
			  # 'font.family': 'Arial',
	          'legend.fontsize': 9,
			  'font.size': 10,
	          'figure.figsize': [6.5,3],
			  'axes.grid': True}
	pylab.rcParams.update(params)			

	#manually read in the data files
	scyllo_15 = numpy.genfromtxt('klv_mon_hR_scyllo_eed_0.2ang_werror.hist')[0:201,:]
	chiro_15 = numpy.genfromtxt('klv_mon_hR_chiro_eed_0.2ang_w_error.hist')[0:201,:]
	scyllo_2 = numpy.genfromtxt('scyllo_eed_withError.histo')[0:201,:]
	chiro_2 = numpy.genfromtxt('chiro_withError_0.2ang.histo1000')[0:201,:]
	water = numpy.genfromtxt('water_eed_withError_0.2ang.histo1000')[0:201,:]

	# combine into two matrices
	high_ratio = numpy.hstack([scyllo_15, chiro_15, water]) # 15:1
	low_ratio = numpy.hstack([scyllo_2, chiro_2, water]) # 2:1

	# create a single figure containing two subplots
	fig = pylab.figure()

	#set the spacing between the two subplots to have no gap between them
	fig.subplots_adjust(top=0.95, bottom=0.15, left=0.1, right=0.95, wspace=0.3)
	ax1 = fig.add_subplot(1,2,1)
	ax2 = fig.add_subplot(1,2,2)

	# turn off the y labels because they share the same y-axis scale
	# yticklabels = ax2.get_yticklabels() + ax2.get_xticklabels()
	yticklabels = ax2.get_xticklabels()
	pylab.setp(yticklabels, visible=False)

	# set a single set of labels for the axes
	ax1.set_xlabel(r'End to end distance ($\AA$)')
	ax1.set_ylabel("Probability")

	isomer_label = ['scyllo', 'chiro', 'water']
	nax=0
	for eed_data in [high_ratio, low_ratio]:
		ax = fig.get_axes()[nax]  #0 = 121, 1 = 122
		ax.set_xlim(0, 23)

		nrows, ncols = eed_data.shape
		for col in range(1, ncols, 3):
			label_index = int(math.ceil((col - 1)/3))
			print "plotting", col, "with error", col+1, label_index
			x = eed_data[:, 0]
			y = eed_data[:, col]
			yerror = eed_data[:, col + 1]
			isomer = isomer_label[label_index]
			ax = plot_with_error(x, y, yerror, config.LABEL[isomer], config.LINE_COLOR[isomer], config.SHADED_COLOR[isomer], axes=ax)
		nax += 1

	ax.legend()
	fig.savefig(file_name, dpi=300)

def main():
	"""docstring for main"""
	plot_eed('KLV_eed_15to1_2to1_with_error.pdf')

if __name__ == '__main__':
	main()


