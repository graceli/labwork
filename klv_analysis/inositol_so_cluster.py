#!/usr/bin/python

import numpy

def average_cluster_columns(datafile, outputname):
    """ this function assumes that every odd column starting from 1 is time or is not used"""
    data = numpy.genfromtxt(datafile)
    sdata = data[0:140001,1::2]
    print sdata
    averaged_columns = numpy.average(sdata, axis=1)
    print averaged_columns
    numpy.savetxt(outputname, numpy.transpose([data[0:140001,0],averaged_columns]), fmt='%0.2f')


average_cluster_columns('scyllo_45to4_nclust_140ns.dat', 'scyllo_45to4_nclust_140ns.txt')
average_cluster_columns('chiro_45to4_nclust_140ns.dat', 'chiro_45to4_nclust_140ns.txt')
