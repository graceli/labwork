import tables
import numpy
import csv

def intersect(sys, nonpolar, polar):
    # get polar matrix
    nrows, ncols = nonpolar.shape

    p_and_np_count  = 0
    nonpolar_count = 0
    polar_count = 0
    total_inositol = 0

    counts = {'sys' : sys, 'polar_only' : 0, 'nonpolar_only' : 0, 'polar_nonpolar' : 0}
    for r in range(0, nrows):
        for i in range(0, ncols):
            if nonpolar[r][i] > 0.0 and polar[r][i] > 0.0:
                counts['polar_nonpolar'] += 1
            elif nonpolar[r][i] > 0.0 and polar[r][i] == 0.0:
                counts['polar_only'] += 1
            elif nonpolar[r][i] == 0.0 and polar[r][i] > 0.0:
                counts['nonpolar_only'] += 1
            total_inositol += 1

    total = counts['polar_nonpolar'] + counts['polar_only'] + counts['nonpolar_only']
    counts['polar_nonpolar'] = counts['polar_nonpolar'] / float(total)
    counts['polar_only'] = counts['polar_only'] / float(total)
    counts['nonpolar_only'] = counts['nonpolar_only'] / float(total)

    return counts


if __name__ == '__main__':
    for ratio in [64]:
        for isomer in ["glycerol", "glucose"]:
            nonpolar_h5file = tables.openFile("nonpolar_contacts_" + str(isomer) + "_" + str(ratio) + ".h5", 'a')
            polar_h5file = tables.openFile("hbonds_" + str(isomer) + "_" + str(ratio) + ".h5", 'a')

            writer = csv.DictWriter(open('intersection_%(ratio)s_%(isomer)s.csv' % vars(), 'wb'), ["sys", "polar_only", "nonpolar_only", "polar_nonpolar"])
            for sys in range(10):
                nonpolar_name = isomer + '_' + str(ratio) + '_' + 'inositol_nonpolar_contacts' + '_' + str(sys)
                polar_name = isomer + '_' + str(ratio) + '_' + 'inositol_hbonds' + '_' + str(sys)
                
                # get nonpolar matrix
                nonpolar = nonpolar_h5file.getNode('/' + nonpolar_name).read().view(dtype=numpy.float64)
                polar = polar_h5file.getNode('/' + polar_name).read().view(dtype=numpy.float64)
                counts = intersect(sys, nonpolar, polar)
                writer.writerow(counts)
