import numpy
import pylab
import glob
import sys
import os

# my module
import utils

def plot_rmsd(iso, ratio, fig, subplot_num):
	#filename="ab_%(iso)s_%(ratio)s_%(subplot_num)s_nosol_whole_rmsd_backbone.xvg" % vars()
	filename="%(iso)s_%(ratio)s_%(subplot_num)s_rmsd.xvg" % vars()
	if os.path.exists(filename):
		data=numpy.genfromtxt(filename)
		nrows,ncols=data.shape
		ax=fig.add_subplot(3, 4, subplot_num)
		fig.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.95, wspace=0.4, hspace=0.4)
		ax.set_xlim(0, 200)
		ax.set_ylim(0, 0.8)
		# ax.axis('tight')
		ax.plot(data[:,0]/1000.0, data[:,1], label=subplot_num)
		ax.plot(data[:,0]/1000.0, [0.5]*nrows, 'g')
		ax.set_title(subplot_num)
		return 0
	else:
		print "WARNING: %(filename)s is not found in the current directory" % vars()
		return 1	

def plot_rmsf(iso, ratio, fig, subplot_num):
	filename="%(iso)s_%(ratio)s_%(subplot_num)d_rmsf.xvg" % vars()
	if os.path.exists(filename):
		print filename
		data=numpy.genfromtxt(filename)
		nrows,ncols=data.shape
		values = data[:,1]
		values_chain_matrix = values.reshape(-1, len(values)/5)
		rows, cols = values_chain_matrix.shape
		# average_rmsf_per_residue = numpy.average(values_chain_matrix, axis=0)
		# pylab.axis('off')

		ax=fig.add_subplot(3, 4, subplot_num)
		fig.subplots_adjust(bottom=0.05,top=0.95,left=0.05,right=0.95,wspace=0.4, hspace=0.4)
		for i in range(0, rows):
			# pylab.plot(values_chain_matrix[i], label=nameslist[fig_num-1]+' chain %(i)d' % vars())
			ax.plot(values_chain_matrix[i], label='chain %(i)d' % vars())
	
		# pylab.xticks(numpy.arange(0,26), matshow_axis_polar())
		# ax.set_yticks()
		ax.grid(True)
		ax.set_xlim(0,26)
		ax.set_ylim(0,0.8)
		ax.set_title(subplot_num)
		# pylab.xlabel('Abeta1-42 residues')
		# pylab.ylabel('RMSF (nm)')
		# pylab.legend(loc=0)
		return 0
	else:
		print "WARNING: %(filename)s is not found in the current directory" % vars()
		return 1

# Some code I wrote to aggregate chain_hbond data
# def comment_xvgr():
# 	# comment out the xvgr header lines for the chain chain hbond analysis
# 	for i in range(1,11):
# 		for chain in range(0,4):
# 			next_chain = chain + 1
# 			os.system("sed -e 's/@/#/g' %(i)d/chain_%(chain)d_%(next_chain)d_hbonds.xvg > chch_hbonds/sys%(i)d_chain_%(chain)d_%(next_chain)d_hbonds.dat" % vars())
# 
# 	# delete the data files for very short (failed) simulation trajs
# 	files = glob.glob('*.dat')
# 	for f in files:
# 		size = os.path.getsize(f)/(1024.0*1024)
# 		if size < 1:
# 			os.system('rm %(f)s' % vars())
# 
# # Example of how its ran
# # for iso in chiro scyllo water; do cd $iso/analysis/hbonds; python ~/Desktop/biophysical_2011/analysis/chain_hbonds.py ../../../${iso}_chain_hbonds_summary.txt; cd ../../../;done
# def chain_hbonds(subplot_num):
# 	average = []
# 	for chain in range(0,4):
# 		next = chain + 1
# 		filename = "%(subplot_num)d/chain_%(chain)d_%(next)d_hbonds.xvg" % vars()
# 
# 		if os.path.exists(filename):
# 			data = numpy.genfromtxt(filename, skip_header=20)
# 			values = data[:,1]
# 			average_hbonds = numpy.average(values)
# 			# print filename, "chain", chain, next, average_hbonds
# 			average.append(average_hbonds)
# 
# 			# ax=fig.add_subplot(3, 4, subplot_num)
# 			# fig.subplots_adjust(bottom=0.05,top=0.95,left=0.05,right=0.95,wspace=0.4, hspace=0.4)
# 			# for i in range(0, rows):
# 			# 	# pylab.plot(values_chain_matrix[i], label=nameslist[fig_num-1]+' chain %(i)d' % vars())
# 			# 	ax.plot(values_chain_matrix[i], label='chain %(i)d' % vars())
# 			# 	
# 			# # pylab.xticks(numpy.arange(0,26), matshow_axis_polar())
# 			# # ax.set_yticks()
# 			# ax.grid(True)
# 			# ax.set_xlim(0,26)
# 			# ax.set_ylim(0,0.8)
# 			# ax.set_title(subplot_num)
# 	print len(average)
# 	if len(average) == 0:
# 		return [0]*4
# 
# 	return numpy.array(average)

# Plots the number of chain-chain hydrogen bonds vs time
def plot_chain_hbond(iso, ratio, fig, subplot_num):
	ax=fig.add_subplot(3, 4, subplot_num)
	fig.subplots_adjust(bottom=0.05,top=0.95,left=0.05,right=0.95,wspace=0.4, hspace=0.4)

	# pylab.rcParams['xtick.labelsize']='10'
	# pylab.rcParams['ytick.labelsize']='12'
	pylab.rcParams['legend.fontsize']='6'
	# pylab.rcParams['figure.figsize'] = [2.5,2.5]
	# pylab.rcParams["axes.titlesize"]='small'

	A=utils.smooth(numpy.genfromtxt('%(iso)s_sys%(subplot_num)d_chain_0_1_hbonds.dat' % vars(), comments="#"), 1000, timestep=2)
	B=utils.smooth(numpy.genfromtxt('%(iso)s_sys%(subplot_num)d_chain_1_2_hbonds.dat' % vars(), comments="#"), 1000, timestep=2) 
	C=utils.smooth(numpy.genfromtxt('%(iso)s_sys%(subplot_num)d_chain_2_3_hbonds.dat' % vars(), comments="#"), 1000, timestep=2)
	D=utils.smooth(numpy.genfromtxt('%(iso)s_sys%(subplot_num)d_chain_3_4_hbonds.dat' % vars(), comments="#"), 1000, timestep=2)

	if A == []:
		print "WARNING: data file empty ... quitting ..."
		return
	
	output_data = numpy.hstack((A, B, C, D))
	numpy.savetxt('%(iso)s_sys%(subplot_num)d_chain_hbonds_1000.dat' % vars(), output_data)
	
	
	time = A[:,0]/1000.0
	ax.plot(time, A[:,1], label="chain 1-2")
	ax.plot(time, B[:,1], label="chain 2-3")
	ax.plot(time, C[:,1], label="chain 3-4")
	ax.plot(time, D[:,1], label="chain 4-5")
	pylab.ylim(0, 30)
	pylab.xlim(0, 200)
	pylab.grid(True)
	# xlabel('Time (ns)')
	# ylabel('Number of interchain hydrogen bonds')
	pylab.legend(loc='lower right')
	pylab.savefig('test_%(iso)s_%(ratio)s_chain_hbond.png' % vars())
		
def main():
	if len(sys.argv) < 3: 
		print "plot_array.py ratio iso analysis"
		sys.exit(0)

	pylab.clf()
	pylab.rcParams['xtick.labelsize']='8'
	pylab.rcParams['ytick.labelsize']='8'
	pylab.rcParams['legend.fontsize']='8'
	pylab.rcParams['figure.figsize'] = [8,6]
	pylab.rcParams["axes.titlesize"]='small'
		
	ratio=sys.argv[1]
	iso = sys.argv[2]
	analysis = sys.argv[3]

	num_plotted = 0
	fig = pylab.figure()
	for i in range(1, 10):
		if analysis == "rmsf":
			print "plotting rmsf"
			num_plotted += plot_rmsf(iso, ratio, fig, i)
		elif analysis == "rmsd":
			print "plotting rmsd"
			num_plotted += plot_rmsd(iso, ratio, fig, i)
		elif analysis == "chain_hbond":
			print "plotting chain_hbond"
			plot_chain_hbond(iso, ratio, fig, i)
		else:
			print "no support for this option"

	# save a rasterized image as a draft b/c for fast viewing  
	# for production figures use pdf/eps
	pylab.savefig(analysis + '.png')

if __name__ == '__main__':
	main()
