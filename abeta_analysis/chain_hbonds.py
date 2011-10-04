import sys
import numpy
import os
import glob

def comment_xvgr():
	# comment out the xvgr header lines for the chain chain hbond analysis
	for i in range(1,11):
		for chain in range(0,4):
			next_chain = chain + 1
			os.system("sed -e 's/@/#/g' %(i)d/chain_%(chain)d_%(next_chain)d_hbonds.xvg > chch_hbonds/sys%(i)d_chain_%(chain)d_%(next_chain)d_hbonds.dat" % vars())

	# delete the data files for very short (failed) simulation trajs
	files = glob.glob('*.dat')
	for f in files:
		size = os.path.getsize(f)/(1024.0*1024)
		if size < 1:
			os.system('rm %(f)s' % vars())

# Example of how its ran
# for iso in chiro scyllo water; do cd $iso/analysis/hbonds; python ~/Desktop/biophysical_2011/analysis/chain_hbonds.py ../../../${iso}_chain_hbonds_summary.txt; cd ../../../;done
def chain_hbonds(subplot_num):
	average = []
	for chain in range(0,4):
		next = chain + 1
		filename = "%(subplot_num)d/chain_%(chain)d_%(next)d_hbonds.xvg" % vars()
	
		if os.path.exists(filename):
			data = numpy.genfromtxt(filename, skip_header=20)
			values = data[:,1]
			average_hbonds = numpy.average(values)
			# print filename, "chain", chain, next, average_hbonds
			average.append(average_hbonds)

			# ax=fig.add_subplot(3, 4, subplot_num)
			# fig.subplots_adjust(bottom=0.05,top=0.95,left=0.05,right=0.95,wspace=0.4, hspace=0.4)
			# for i in range(0, rows):
			# 	# pylab.plot(values_chain_matrix[i], label=nameslist[fig_num-1]+' chain %(i)d' % vars())
			# 	ax.plot(values_chain_matrix[i], label='chain %(i)d' % vars())
			# 	
			# # pylab.xticks(numpy.arange(0,26), matshow_axis_polar())
			# # ax.set_yticks()
			# ax.grid(True)
			# ax.set_xlim(0,26)
			# ax.set_ylim(0,0.8)
			# ax.set_title(subplot_num)
	print len(average)
	if len(average) == 0:
		return [0]*4
		
	return numpy.array(average)
	
if __name__ == '__main__':
	
	filename = sys.argv[1]
	data = []
	for i in range(1,11):
		data.append(chain_hbonds(i))
	
	print data
	matrix = numpy.vstack(data)
	print matrix
	numpy.savetxt(filename, matrix, fmt='%.3f')