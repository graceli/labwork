#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
from xvg2png import xvg2array

def define_row_col(len_infiles):
    if len_infiles <= 2 :
        row =2; col = 1
    elif len_infiles <= 4:
        row = col = 2
    elif len_infiles <= 6:
        row = 2; col = 3
    elif len_infiles <= 9:
        row = col = 3
    elif len_infiles <= 12:
        row = 3; col = 4
    elif len_infiles <= 16:
        row = col = 4
    elif len_infiles <= 20:
        row = 4; col = 5
    elif len_infiles <= 25:
        row = col = 5
    elif len_infiles <= 30:
        row = 5; col = 6
    else:
        row = col = None
        print "too many graphs!"
    return row, col

def quickplot_nolegend(infile,ax):
    x, y = xvg2array(infile)
    ax.plot(x,y)
    ax.set_title(infile)
    ax.grid()

def quicksubplots(infiles):
    fig = plt.figure()
    len_infiles = len(infiles)
    row, col = define_row_col(len_infiles)
    ks = range(len(infiles))
    for k, infile in zip(ks, infiles):
        ax = fig.add_subplot(row,col,k+1)
        quickplot_nolegend(infile,ax)

if __name__ == '__main__':
    infiles = sorted(glob.glob('*%s*' % sys.argv[1]))
    quicksubplots(infiles)
    plt.show()

