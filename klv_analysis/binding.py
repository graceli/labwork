import re
import os
import sys
import csv
import glob
import math
import numpy
import tables
import plot_and_save2hdf5 as myh5
import utils

# First look for the code I've written before ... I definitely remember doing this before
# def binding_constant():
    # Add up polar and nonpolar matrices and sum over axis=1
    # Count the number of nonzero elements

    
def intersection_mon(h5file, csv_file, isomer, ratio):
    polar_matrix = myh5.getTableAsMatrix(h5file, '/inositol/inos_total')
    nonpolar_matrix = myh5.getTableAsMatrix(h5file, '/residue/per_inos_contacts')
    
    print polar_matrix.shape
    print nonpolar_matrix.shape
    
    assert polar_matrix.shape == nonpolar_matrix.shape, "the two matrices are expected to have the same dimensions"
    
    nrows, ncols = polar_matrix.shape
    counts = [{'polar_only':0, 'nonpolar_only':0, 'polar_nonpolar':0}, {'polar_only':0, 'nonpolar_only':0, 'polar_nonpolar':0}]
    for i in range(0, nrows): 
        for j in range(1, ncols):
            if polar_matrix[i][0] < 7500:
                if polar_matrix[i][j] and nonpolar_matrix[i][j]:
                    counts[0]['polar_nonpolar'] += 1
                elif polar_matrix[i][j]:
                    counts[0]['polar_only'] += 1
                elif nonpolar_matrix[i][j]:
                    counts[0]['nonpolar_only'] += 1
            else:
                if polar_matrix[i][j] and nonpolar_matrix[i][j]:
                    counts[1]['polar_nonpolar'] += 1
                elif polar_matrix[i][j]:
                    counts[1]['polar_only'] += 1
                elif nonpolar_matrix[i][j]:
                    counts[1]['nonpolar_only'] += 1

    # class csv.DictWriter(csvfile, fieldnames[, restval=''[, extrasaction='raise'[, dialect='excel'[, *args, **kwds]]]])
    total = counts[0]['polar_only'] + counts[0]['nonpolar_only'] + counts[0]['polar_nonpolar']
    fraction = {'polar_only': float(counts[0]['polar_only']) / total, 'nonpolar_only' : float(counts[0]['nonpolar_only']) / total, 'polar_nonpolar' : float(counts[0]['polar_nonpolar']) / total}

    print counts[0]

    writer = csv.DictWriter(open(csv_file, 'wb'), counts[0].keys())
    writer.writeheader()
    writer.writerow(counts[0])
    writer.writerow(fraction)
 

def intersection_disordered(h5file, ratio, system_indices):
    """intersection analysis for disordered oligomers"""

    #nasty fix for different table names
    polarName = {'15to4' : 'whole_nosol_0-200ns_inos_total.dat', '45to4' : 'whole_nosol_0-200_inos_total.dat'}
    nonpolarName = {'15to4' : 'whole_nosol_0-200ns_per_inositol_contacts.dat', '45to4' : 'whole_nosol_0-200_per_inositol_contacts.dat'}
    
    isomerList = ["scyllo", "chiro"]
    polarPath = "/polar"
    nonpolarPath = "/nonpolar_residue"
    dataList = [['isomer', 'system#', 'polar_only', 'polar_and_nonpolar', 'nonpolar_only', 'total']]
    
    resultsWriter = csvwriter.writer(open(ratio + '_intersection.txt', 'wb'), delimiter=' ')
    
    print system_indices
    
    for iso in isomerList:
        data = []
        for sys in system_indices[iso]:
            polarFile = os.path.join(polarPath, "%(iso)s_sys%(sys)s_%(ratio)s_" % vars() + polarName[ratio])
            print "analyzing", polarFile
            polarMatrix = myh5.getTableAsMatrix(h5file, polarFile, dtype=numpy.float64)
            print polarMatrix
            nonpolarFile = os.path.join(nonpolarPath, "%(iso)s_sys%(sys)s_%(ratio)s_" % vars() + nonpolarName[ratio])
            nonpolarMatrix = myh5.getTableAsMatrix(h5file, nonpolarFile, dtype=numpy.float64)[::2, :]
            
            if polarMatrix != None and nonpolarMatrix != None:
                rows, cols = nonpolarMatrix.shape
                print rows, cols
                print polarMatrix.shape
                polar_and_nonpolar = 0.0
                polar_only = 0.0
                nonpolar_only = 0.0
                rows = min(rows, 100000)
                print rows, cols
                for i in range(20001, rows):
                    for j in range(1, cols):
                        if polarMatrix[i][j] and nonpolarMatrix[i][j]:
                            polar_and_nonpolar += 1
                        elif polarMatrix[i][j]:
                            polar_only += 1
                        elif nonpolarMatrix[i][j]:
                            nonpolar_only += 1

                total = polar_only + polar_and_nonpolar + nonpolar_only 
                dataList.append([iso, sys, polar_only / total, polar_and_nonpolar / total, nonpolar_only / total, total])
                data.append([polar_only / total, polar_and_nonpolar / total, nonpolar_only / total, total])
        
        print data
        print numpy.array(data)
        average = numpy.average(numpy.array(data), axis=0)
        std = numpy.std(numpy.array(data), axis=0)
        print average.tolist()
        print std.tolist()
        # 
        addToListAvg = [iso+' avg', 'all']
        addToListAvg.extend(average.tolist())
        addToListStd = [iso+' std', 'all']
        addToListStd.extend(std.tolist())
        
        dataList.append(addToListAvg)
        dataList.append(addToListStd)

    # numpy.savetxt("15to4_intersection.gz", dataList, fmt='%s %d %0.3f %0.3f %0.3f %d')
    resultsWriter.writerows(dataList)

def intersection_beta():
    isomerList = ["scyllo", "chiro"]
    polar_path = "/polar"
    nonpolar_path = "/nonpolar_residue"
    
    for iso in isomerList:
        for sys in range(0, 6):
            nonpolar_file = os.path.join(nonpolar_path, "%(iso)s_t%(sys)d_per_inositol_contacts.dat" % vars())
            polar_file = os.path.join(polar_path, "%(iso)s_t%(sys)d_inos_total.dat" % vars())

            nonpolar_matrix = numpy.genfromtxt(nonpolar_file)
            polar_matrix = numpy.genfromtxt(polar_file)

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

            writer = csv.DictWriter(open('%(iso)s_intersection%(sys)d.csv' % vars(), 'wb'), counts.keys())
            writer.writerow(counts)
            

# This function used to be a script that I used to calculate the nonpolar_residue contact binding bar plots 
# for the monomeric system at a high concentration.
# Note that the concentration is mostly irrelevant for the coding.  
# This was originally written as a script only for computing this subset of monomer data.
# The difference lies in the nature of the dataset. I should extract out the common logic which is that 
# for N system, you calculate a average and a standard deviation.
def nonpolar_residue_monomer_15to1():
    # open X number of files for the monomer (there are 500)
    # chiro_sys249_per_res_contacts.dat
    fileList = glob.glob("*per_res_contacts.dat")
    TOTAL_FILES = len(fileList)
    NRES = 7
    segments = [[0]*NRES, [0]*NRES, [0]*NRES, [0]*NRES, [0]*NRES]
    average = [0]*NRES

    seg = 0
    filesProcessed = 0
    for file in fileList:
        print "processing", file

        r = csv.reader(open(file), delimiter=' ')   
        nrows = 0
        for row in r:
            time = float(row[0])
            for i in range(1,NRES+1):
                value = int(row[i])
                segments[seg][i-1] += value
                average[i-1] += value
            nrows += 1

        filesProcessed += 1
        print "processed", filesProcessed

        if filesProcessed == 100:
            seg += 1    
            print "at seg", seg
            filesProcessed = 0

    # compute the standard deviation for each value in the array 
    sd = [0]*7
    TOTAL_FRAMES = nrows
    NSEGS = seg

    print "total frames", TOTAL_FRAMES, "NSEGS", seg

    for i in range(0, NRES):
        for s in range(0, seg):
            sd[i] += math.pow(float(segments[s][i])/(100*TOTAL_FRAMES) - float(average[i])/(TOTAL_FILES*TOTAL_FRAMES),2)
        sd[i] = math.sqrt(sd[i]/NSEGS)

    # output the file 
    # one column for the final values, second column is the stdev
    print "# residue average sd"
    print "# files Processed", filesProcessed
    print "# total frames per file", TOTAL_FRAMES
    print "$ total segments", seg
    for i in range(0,NRES):
        print i, float(average[i]) / (TOTAL_FILES * TOTAL_FRAMES), sd[i]


def nonpolar_residue_disordered(h5file, tag):
    scyllo_pattern = re.compile(r'scyllo')
    chiro_pattern = re.compile(r'chiro')
    atype_pattern = re.compile(r'residue_contact')
    data_list = {'scyllo':[], 'chiro':[]}

    #fix this number for now
    N_datapoints = 190000
    for table in h5file.listNodes("/nonpolar_residue", 'Table'):
        if atype_pattern.search(table.name):
            table_path = os.path.join("/nonpolar_residue", table.name)
            data = myh5.getTableAsMatrix(h5file, table_path, dtype=numpy.float64)
            
            print data
            
            sum_over_time = numpy.average(data[20000:N_datapoints, 1:], axis = 0)
            
            print sum_over_time
            
            # This matrix is Nres by 4, where 4 is the number of peptides in the system (disordered oligomer)
            # Each row of this matrix represents a single amino acid
            # Each column is a peptide sequence A, E, L, K, F, F, V
            sum_over_time.shape = (sum_over_time.size / 4, 4)
            
            # Average over all peptides in the system (over columns, hence axis = 1). 
            # The resulting array of numbers has units of per peptide.
            avg_over_peptides = numpy.average(sum_over_time, axis = 1)
            
            if scyllo_pattern.search(table.name):
                data_list['scyllo'].append(avg_over_peptides)
            elif chiro_pattern.search(table.name):
                data_list['chiro'].append(avg_over_peptides)
            else:
                print "No pattern matches", table.name

    # save results to flat files
    for isomer in data_list.keys():
        nparray = numpy.array(data_list[isomer])
            
        # dump the list of counts for each system
        numpy.savetxt('%(tag)s_nonpolar_residue_inositol_contact_%(isomer)s_counts.txt' % vars(), nparray, fmt='%0.8f')
        print "saved", isomer, "analysis with shape", nparray.shape

        # average over all the systems; each system is a row in nparray
        average = numpy.average(nparray, axis=0)
        std = numpy.std(nparray, axis=0)/math.sqrt(8)
        numpy.savetxt('%(tag)s_nonpolar_residue_inositol_contact_%(isomer)s_avg_std.txt' % vars(), [average, std], fmt='%0.8f')

# This was originally the script binding_pattern.py which operates on the beta oligomer nonpolar residue
# binding dataset
def nonpolar_residue_beta(isomer):
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

if __name__ == '__main__':
    if len(sys.argv) < 3:
        parser.error("Please specify a .h5 input file and a tag for output")
    
    filename = sys.argv[1]
    tag = sys.argv[2]
    # option = sys.argv[2]
    # use_flat_flag = False
    print filename
    
    h5file = tables.openFile(filename)
    # isomer,sys,analysis = filename.split('_')
    nonpolar_residue_disordered(h5file, tag)
    # intersection_mon(h5file, isomer, '15to1') 
