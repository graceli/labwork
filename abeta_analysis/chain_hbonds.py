import sys
import numpy
import os

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