import numpy
import tables
import pylab

def plot_contact_map(data, filename):
    data = numpy.genfromtxt(filename)
    pylab.matshow(data)
    pylab.clim(limits)
    pylab.colorbar(orientation="horizontal", shrink=1)

    # print "the plot axis labels:", data_axis  

    #how the axis will be rendered
    x = numpy.arange(0, len(data[0]))
    xlabels = tuple(data_axis)
    pylab.xticks(x, xlabels)

    y = numpy.arange(0, len(data[:,0]))
    ylabels = tuple( ["chain 1", "chain 2", "chain 3", "chain 4", "chain 5"] )
    pylab.yticks(y, ylabels)

    #save the figure to disk as a pdf
    pylab.savefig(filename + '.pdf' % vars())   
    numpy.savetxt(filename + '.dat', numpy.transpose(data), fmt='%f')



polar_data_axis = [ 'L', 'V', 'F', 'F', 'A', 'E', 'D', 'V', 'G', 'S', 'N', 'K', 'G', 'A', 'I', 'I', 'G', 'L', 'M', 'V', 'G', 'G', 'V', 'V', 'I', 'A' ]
nonpolar_data_axis = [ 'L', 'V', 'F', 'F', 'A', 'E', 'D', 'V', 'S', 'N', 'K', 'A', 'I', 'I', 'L', 'M', 'V', 'V', 'V', 'I', 'A' ]

for ratio in [15, 64]:
    for isomer in [ "scyllo", "chiro", "glycerol", "glucose" ]:
        # nonpolar_data = tables.openFile(str(ratio) + "_" + isomer + ".h5" ).getNode('/' + ...)
        # hbonds_data = tables.openFile(str(ratio) + "_" + isomer + ".h5" ).getNode('/' + ...)

        nonpolar_name = isomer + "_" + str(ratio) + "_" + "nonpolar_contact_matrix.txt"
        hbonds_name = isomer + "_" + str(ratio) + "_" + "hbonds_contact_matrix.txt"
        nonpolar_data = numpy.genfromtxt(nonpolar_name)
        hbonds_data = numpy.genfromtxt(hbonds_name) 
        plot_contact_map(nonpolar_data, nonpolar_name)
        plot_contact_map(hbonds_data, hbonds_name)
