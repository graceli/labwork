import numpy
import pylab
import os
import math
import config
import config_plot

def plot_rmsd_with_std_err(iso, ratio):
    total = []
    count = 0
    # Average data over all replicas
    length = 110000
    if ratio == 64 and (iso == "chiro" or iso == "scyllo"):
        length = 70000
    if iso == "glucose":
        length = 100000
    for subplot_num in range(10):     
        filename="%(ratio)s/%(iso)s/rmsd/%(subplot_num)s_rmsd.xvg" % vars()
        print filename
        if os.path.exists(filename):
            data = numpy.genfromtxt(filename)[0:length]
            print "appended data with dimensions:", data.shape
            total.append(data)
            count += 1
        else:
            print filename, "does not exist"
            
    nd_total = numpy.dstack(total)
    print nd_total
    print nd_total.shape
    mean_data = numpy.average(nd_total, axis=2)
    print mean_data.shape
    std_error = numpy.std(nd_total, axis=2)/math.sqrt(count - 1)
    print std_error.shape

    fig = pylab.figure()
    ax=fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.95, wspace=0.4, hspace=0.4)
    ax.set_ylim(0, 0.6)
    ax.set_xlim(20, 200)
    x = mean_data[:,0]/1000.0
    y = mean_data[:,1]
    yerror = std_error[:,1]
    ax.fill_between(x, y-yerror, y+yerror, alpha=0.5, facecolor='pink', edgecolor='none')
    ax.plot(x, y)
    ax.plot(x, [0.5]*x.shape[0], 'g')
    fig.savefig(iso + '_' + str(ratio) + 'rmsd_with_stderr.png')

def plot_rmsf(iso, ratio):
    data_all = []    
    for subplot_num in range(10):    
        filename="%(ratio)s/%(iso)s/rmsf_calpha/%(iso)s_%(ratio)s_%(subplot_num)d_rmsf.xvg" % vars()
        if os.path.exists(filename):
            print filename
            data=numpy.genfromtxt(filename)
            nrows,ncols=data.shape
            print data.shape
            values = data[:,1]
            print values.shape
            values_chain_matrix = values.reshape(-1, len(values)/5)
            print values_chain_matrix.shape
            # rows, cols = values_chain_matrix.shape
            data_all.append(values_chain_matrix)
    
    nd_data_all = numpy.dstack(data_all)
    print nd_data_all.shape
    
    average_rmsf = numpy.average(nd_data_all, axis=2)
    print average_rmsf.shape
    rows, cols = average_rmsf.shape
    
    fig = pylab.figure()
    ax=fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.1,top=0.90,left=0.15,right=0.90,wspace=0.0, hspace=0.0)
    
    for i in range(0, rows):
        # pylab.plot(values_chain_matrix[i], label=nameslist[fig_num-1]+' chain %(i)d' % vars())
        ax.plot(average_rmsf[i], label='chain %(i)d' % vars())

    # pylab.xticks(numpy.arange(0,26), matshow_axis_polar())
    # ax.set_yticks()
    # ax.grid(True)
    pylab.axis('tight')
    ax.set_xlim(0,26)
    ax.set_ylim(0,0.7)
    
    pylab.xticks(numpy.arange(0,26), matshow_axis_polar())
    pylab.yticks()
	
    fig.savefig(iso + '_' + str(ratio) + '_avg_rmsf.png')
    # ax.set_title(subplot_num)
    # pylab.xlabel('Abeta1-42 residues')
    # pylab.ylabel('RMSF (nm)')
    # pylab.legend(loc=0)


def plot_all_rmsd():
    for ratio in config.ratio_list:
        for isomer in config.isomer_list:
            print isomer, ratio, "RMSD with standard error"
            plot_rmsd_with_std_err(isomer, ratio)

    # plot_rmsd_with_std_err("water", 15)
    # plot_rmsd_with_std_err("water", 64)
    # plot_rmsd_with_std_err("glucose", 64)
                   
def plot_all_rmsf():
    for ratio in config.ratio_list:
        for isomer in config.isomer_list:
            print isomer, ratio, "plotting average rmsf ..."
            plot_rmsf(isomer, ratio)


def matshow_axis_polar():
	""" hardcoded residue names for the polar data axis (all the residues in Abeta17-42) """

	residue_map = ['L', 'V', 'F', 'F', 'A', 'E', 'D', 'V', 'G', 'S', 'N', 'K',
	'G', 'A', 'I', 'I', 'G', 'L', 'M', 'V','G', 'G', 'V','V', 'I', 'A']

	return residue_map
            	            
def main():
    config_plot.configure_plot()
    
    print "Plotting"
    # plot_all_rmsd()
    # plot_all_rmsf()
    
    config_plot.configure_large_plot()
    plot_rmsf("water", 64)
    

    
if __name__ == '__main__':
    main()