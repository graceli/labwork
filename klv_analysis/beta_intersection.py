import numpy
import csv

def intersection():
    for sys in range(0, 6):
        nonpolar_file = "nonpolar/t%(sys)d_per_inositol_contacts.dat" % vars()
        polar_file = "polar/data/t%(sys)d_inos_total.dat" % vars()

        nonpolar_matrix_full = numpy.genfromtxt(nonpolar_file)
        polar_matrix = numpy.genfromtxt(polar_file)
        nonpolar_matrix = nonpolar_matrix_full[::2,:]

        print polar_matrix.shape
        print nonpolar_matrix.shape

        assert polar_matrix.shape == nonpolar_matrix.shape, "the two matrices are expected to have the same dimensions"
        
        nrows, ncols = polar_matrix.shape
        counts = {'polar_only' : 0, 'nonpolar_only' : 0, 'polar_nonpolar' : 0}
        total = 0.0
        for i in range(1, nrows):
            for j in range(1, ncols):
                total += 1.0
                if polar_matrix[i][j] and nonpolar_matrix[i][j]:
                    counts['polar_nonpolar'] += 1
                elif polar_matrix[i][j]:
                    counts['polar_only'] += 1
                elif nonpolar_matrix[i][j]:
                    counts['nonpolar_only'] += 1

        # normalize
        total = counts['polar_nonpolar'] + counts['polar_only'] + counts['nonpolar_only']
        counts['polar_nonpolar'] = counts['polar_nonpolar'] / float(total)
        counts['polar_only'] = counts['polar_only'] / float(total)
        counts['nonpolar_only'] = counts['nonpolar_only'] / float(total)
        
        # class csv.DictWriter(csvfile, fieldnames[, restval=''[, extrasaction='raise'[, dialect='excel'[, *args, **kwds]]]])
        print sys, counts

        writer = csv.DictWriter(open('intersection%(sys)d.csv' % vars(), 'wb'), counts.keys())
        writer.writerow(counts)


intersection()
