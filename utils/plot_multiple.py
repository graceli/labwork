import sys
import os
import glob
import pylab
import math
import numpy

def plot(filename):
	"""docstring for plot"""
	data_array = numpy.genfromtxt(filename)
	pylab.plot(data_array[:,0], data_array[:,1], linewidth=2.0)
	
	# line, = plt.plot(x, y, '-')
	# line.set_antialiased(False) # turn off antialising
	
def main():
	"""docstring for main"""
	bash_expr = sys.argv[1]
	files = glob.glob(bash_expr)
	num_files = len(files)
	print "You are plotting", num_files
	if num_files > 25:
		print "You've requested to plot", num_files, "Are you sure?"
	
	cols = math.ceil(math.sqrt(num_files))
	rows = math.ceil(num_files/float(cols))
	
	print "This will generate a plot with",rows, "rows and", cols, "cols"
	
	counter=1
	for r in range(1, int(rows)+1):
		for c in range(1, int(cols)+1):
			file_to_plot = files[counter-1]
			print "plot", file_to_plot, "in row", r, "and col", c, "number", counter
			pylab.subplot(int(rows), int(cols), counter)
			plot(file_to_plot)
			counter += 1
			if counter > num_files:
				break
				
	pylab.show()
	
if __name__ == '__main__':
	main()