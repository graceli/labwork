import subprocess
import os
import numpy
import tables
import tarfile
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


def process_hbonds_inositol(h5file):
    for isomer in config.isomer_list:
        for ratio in config.ratio_list:
            tar = tarfile.open(os.path.join(config.data_source, "analysis_{0}_{1}_hbonds_inositol.tgz".format(isomer, ratio)))
            
            for sys in range(10):
                data_path = os.path.join('hbonds_inositol', str(sys))

                # Each file here represents a single inositol molecule and 
                # the number of hbonds it froms with the aggregate
                # Form a matrix with column ordering of 1 ... 64 and save this
                files = [ os.path.join(data_path, str(i) + '.xvg') for i in range(1, ratio+1) ]
                first_file = True
                column_stack = []
                for f in files:
                    print "Extracting and reading file", f
                    data = numpy.genfromtxt(tar.extractfile(f))
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
                table_name = isomer + '_' + str(ratio) + '_' + 'hbonds_inositol' + '_' + str(sys)
                # Each item in list is a row vector, stack them vertically and transpose into a matrix with dimensions time X Ninositols 
                data_all = numpy.transpose(numpy.vstack(column_stack))

                print "Saving data from", files, "into", h5file.root, table_name
                atom = tables.Atom.from_dtype(data_all.dtype)
                filters = tables.Filters(complib='zlib', complevel = 5)
                h5_carray = h5file.createCArray(h5file.root, table_name, atom, data_all.shape, filters=filters)
                h5_carray[:] = data_all[10000:,:]

def process_nonpolar_inositol(h5file):
    """docstring for process_nonpolar_inositol"""
    # read in each of the chains and sum up the matrices
    for isomer in config.isomer_list:
        for ratio in config.ratio_list:
            tar = tarfile.open(os.path.join(config.data_source, "analysis_{0}_{1}_nonpolar.tgz".format(isomer, ratio)))
            for sys in range(10):
                # Sum data files from all chains
                # Each file contains the inositol nonpolar contacts to a single chain in the protofibril
                # ie. for each inositol molecule is it bound to chain X, where X=1,2,3,4,5?
                # Note that the sum of all of these matrices is the total number of 
                # nonpolar contacts made by inositol to the protofibril
                data_all = None
                for ch in range(5):
                    data_file = os.path.join('nonpolar', '%(sys)s_chain%(ch)s_inositol_np_contact' % vars() + '.dat')
                    # extract the file from the tar as a file object
                    member = tar.extractfile(data_file)
                    a_chain_data = numpy.genfromtxt(member)

                    if ch == 0:
                        data_all = a_chain_data
                    else:
                        data_all[:, 1:] += a_chain_data[:,1:]

                table_name = isomer + '_' + str(ratio) + '_' + 'nonpolar' + '_' + str(sys)

                print "Saving data into", h5file.root, table_name

                # Saves into h5 file using pytables as a CArray object. Code straight out of documentation
                # http://readthedocs.org/docs/pytables/en/latest/usersguide/libref.html?highlight=carray#tables.CArray
                # CArray chosen for homogeneous data
                atom = tables.Atom.from_dtype(data_all.dtype)
                filters = tables.Filters(complib='zlib', complevel = 5)
                h5_carray = h5file.createCArray(h5file.root, table_name, atom, data_all.shape, filters=filters)
                h5_carray[:] = data_all[9900:,:]
                
def main():
    h5file = tables.openFile("analysis.h5", 'w')
    print "reading"
    process_hbonds_inositol(h5file)
    process_nonpolar_inositol(h5file)
    h5file.close()

if __name__ == '__main__':
    main()
