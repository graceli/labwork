import numpy
import math

def nonpolar_residue(isomer):
    data_list = []

    for i in range(0, 6):
        file = "%(isomer)s_t%(i)d_per_residue_contacts.dat" % vars()
        print "analyzing ", file

        # read in the file
        data = numpy.genfromtxt(file, comments="#", dtype='float')
        nrows,ncols = data.shape
        print nrows

        # sum over rows
        time_avg = numpy.average(data[:,1:], axis=0)
        print time_avg.shape
        time_avg.shape = (time_avg.size / 16, 16) 
        print time_avg.shape
        sum_over_peptides = numpy.sum(time_avg, axis = 1)
                
        data_list.append(sum_over_peptides)

    # save results to flat files
    nparray = numpy.array(data_list)
                
    # dump the list of results for each system
    numpy.savetxt('%(isomer)s_nonpolar_residue_contact.txt' % vars(), nparray, fmt='%0.8f')

    # average over all the systems; each system is a row in nparray
    average = numpy.average(nparray, axis=0) / 16 
    std = numpy.std(nparray, axis=0) / 16 / math.sqrt(6)

    #save the normalized average and std
    numpy.savetxt('%(isomer)s_nonpolar_residue_contact_avg_std.txt' % vars(), [average, std], fmt='%0.8f')

nonpolar_residue("scyllo")
nonpolar_residue("chiro")
 
