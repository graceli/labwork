import subprocess
import os
import tarfile
import sys
import numpy
import tables
import config

# Note: all my operations involves the following
# TODO can I build a system (or is there an existing one) to suit my needs?
# Analysis -> tar.gz files with a certain structure
# Data munging app -> read tar.gz files with these known structures -> performs custom data integration into a database (pytables) for downstream analysis
# Analysis app -> performs calculations on the munged data, either outputs a csv file with results and or a plot corresponding to the results

# This script solely works with tarfiles and entirely eliminates working with flat files on disk
# I think this is super efficient
# 1) Disk space efficiency
# 2) Reduces the number of flat ascii files
# 3) Speeds up processing by touching disk less (don't have to extract files to disk, then read from disk)
#     Could even copy files to memory and then process files from /dev/shm


def process_hbonds(h5file, isomer, ratio, sys_idx, analysis_type, num_residues):
    print "process_hbonds: munging files for", analysis_type, "with", num_residues
    tar_file_name = "analysis_{0}_{1}_hbonds_inositol.tgz".format(isomer, ratio)
    if analysis_type == "residue":
        tar_file_name = "analysis_{0}_{1}_hbonds.tgz".format(isomer, ratio)
    
    tar_file_path = os.path.join(config.data_source, tar_file_name)
    
    if not os.path.exists(tar_file_path):
        print "copying", tar_file_name, "to", config.data_source
        os.system("cp {0} {1}".format(tar_file_name, config.data_source))

    tar = tarfile.open(tar_file_path)
    
    data_path = os.path.join('hbonds_inositol', str(sys_idx))
    if analysis_type == "residue":
        data_path = os.path.join('hbonds', str(sys_idx))
    
    # Each file here represents a single inositol molecule and 
    # the number of hbonds it froms with the aggregate
    # Form a matrix with column ordering of 1 ... 64 and save this
    files = [os.path.join(data_path, str(i) + '.xvg') for i in range(1, num_residues+1)]
    if analysis_type == "residue":
        files = [os.path.join(data_path, str(i) + '.xvg') for i in range(0, num_residues)]
    
    first_file = True
    column_stack = []
    for f in files:
        print "Extracting and reading file", f, "from", tar_file_path
        data = numpy.genfromtxt(tar.extractfile(f), comments="#")
        print data.shape

        if first_file:
            # grab the first two columns
            time = data[:, 0]
            column_stack.append(time)
            column_stack.append(data[:,1])
            first_file = False
        else:
            # only keep the 1st column (get rid of 0th, and 2nd column)
            column_stack.append(data[:, 1])

    # convert the list into a numpy matrix and save into pytables
    table_name = isomer + '_' + str(ratio) + '_' + 'hbonds_inositol' + '_' + str(sys_idx)
    # Each item in list is a row vector, stack them vertically and transpose into a matrix with dimensions time X Ninositols 
    data_all = numpy.transpose(numpy.vstack(column_stack))
    
    print "Saving data from", files, "into", h5file.root, table_name
    atom = tables.Atom.from_dtype(data_all.dtype)
    filters = tables.Filters(complib='zlib', complevel = 5)
    h5_carray = h5file.createCArray(h5file.root, table_name, atom, data_all.shape, filters=filters)
    h5_carray[:] = data_all

def process_nonpolar(h5file, isomer, ratio, sys_idx, analysis_type):
    print "process_nonpolar: munging files for", analysis_type

    # read in each of the chains and sum up the matrices
    tar = tarfile.open(os.path.join(config.data_source, "analysis_{0}_{1}_nonpolar.tgz".format(isomer, ratio)))
    # Sum data files from all chains
    # Each file contains the inositol nonpolar contacts to a single chain in the protofibril
    # ie. for each inositol molecule is it bound to chain X, where X=1,2,3,4,5?
    # Note that the sum of all of these matrices is the total number of 
    # nonpolar contacts made by inositol to the protofibril
    data_all = None
    for ch in range(5):
        data_file = os.path.join('nonpolar', '%(sys_idx)s_chain%(ch)s_%(analysis_type)s_np_contact' % vars() + '.dat')
        
        # extract the file from the tar as a file object
        print "Extracting", data_file
        
        member = tar.extractfile(data_file)
        a_chain_data = numpy.genfromtxt(member, comments="#")

        if ch == 0:
            data_all = a_chain_data
        else:
            data_all[:, 1:] += a_chain_data[:,1:]

    table_name = isomer + '_' + str(ratio) + '_' + 'nonpolar' + '_' + str(sys_idx)

    print "Saving data into", h5file.root, table_name

    # Saves into h5 file using pytables as a CArray object. Code straight out of documentation
    # http://readthedocs.org/docs/pytables/en/latest/usersguide/libref.html?highlight=carray#tables.CArray
    # CArray chosen for homogeneous data
    atom = tables.Atom.from_dtype(data_all.dtype)
    filters = tables.Filters(complib='zlib', complevel = 5)
    h5_carray = h5file.createCArray(h5file.root, table_name, atom, data_all.shape, filters=filters)
    h5_carray[:] = data_all
                
def main():
    if len(sys.argv) < 2:
        print "usage: clean.py isomer ratio sys"
        sys.exit(0)

    isomer = sys.argv[1]
    ratio = int(sys.argv[2])
    sys_idx = int(sys.argv[3])
 
    h5file = tables.openFile("analysis.h5", 'w')

    print "Munging analysis files into pytables"
    
    # reading hbonds for all inositol molecules
    process_hbonds(h5file, isomer, ratio, sys_idx,  "inositol", 15)
    
    # reading hbonds for all protein residues
    process_hbonds(h5file, isomer, ratio, sys_idx, "residue", 130)
    
    # nonpolar contact matrix for inositol
    process_nonpolar(h5file, isomer, ratio, sys_idx, "inositol")
    
    # nonpolar contact matrix for protein residues
    process_nonpolar(h5file, isomer, ratio, sys_idx, "residue")
    
    h5file.close()

if __name__ == '__main__':
    main()
