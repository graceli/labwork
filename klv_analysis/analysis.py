#!/usr/bin/env python
import os
import sys
import re
import csv
import glob
import subprocess
from optparse import OptionParser

import pylab
import numpy
import tables

#how to share imports?
import plot_and_save2hdf5 as myh5
import utils

# my analysis module code
# import dssp
# import binding
import timeseries_hb
import timeseries_nonpolar
import cluster
import config

# implement using the Strategy Pattern?


# inositol_so_cluster.py
# def average_cluster_columns(datafile, outputname):
#     """ this function assumes that every odd column starting from 1 is time or is not used"""
#     data = numpy.genfromtxt(datafile)
#     sdata = data[0:140001,1::2]
#     print sdata
#     averaged_columns = numpy.average(sdata, axis=1)
#     print averaged_columns
#     numpy.savetxt(outputname, numpy.transpose([data[0:140001,0],averaged_columns]), fmt='%0.2f')

# average_cluster_columns('scyllo_45to4_nclust_140ns.dat', 'scyllo_45to4_nclust_140ns.txt')
# average_cluster_columns('chiro_45to4_nclust_140ns.dat', 'chiro_45to4_nclust_140ns.txt')

def main():
	parser = OptionParser()
	parser.add_option("-f", "--use-flat", action="store_true", dest="use_flat_flag", default=False,
						help="load plot data from flat files")
	parser.add_option("-n", "--nonpolar-timeseries", action="store_true", dest="run_nonpolar_flag", default=False,
	                  help="run nonpolar timeseries analysis and plot")
	parser.add_option("-p", "--hb-timeseries", action="store_true",
	                  dest="run_polar_flag", default=False,
	                  help="run polar timeseries analysis and plot")
	parser.add_option("-c", "--cluster", action="store_true", dest="run_cluster", default=False,
						help="process and plot cluster size analysis")

	(options, args) = parser.parse_args()

	# print options
	# print args
	if len(args) < 1:
		parser.error("Please specify a .h5 input file")
	
	filename = args[0]
	# option = sys.argv[2]
	# use_flat_flag = False

	h5file = tables.openFile(filename)
	ratio,ext = os.path.splitext(filename)

	# binding.nonpolar_residue(h5file, ratio)
	# binding.intersection(h5file, ratio)
	# dssp.run(h5file)
	
	config.configure_plot()

	if options.run_polar_flag:
		timeseries_hb.run(h5file, ratio, use_flat_files=options.use_flat_flag)
	
	if options.run_nonpolar_flag:
		timeseries_nonpolar.run(h5file, ratio, use_flat_files=options.use_flat_flag)
		
	if options.run_cluster:
		cluster.run(h5file, ratio, use_flat_files=options.use_flat_flag)

if __name__ == '__main__':
	main()