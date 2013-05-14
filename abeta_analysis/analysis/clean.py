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

class DataMatrix:
    def __init__(self):
        column_stack = []
    
    def append(array):
        # Appends a column to the data matrix
        column_stack.append(array)

    def save_to_h5(h5file, table_name):
        # Each item in list is a row vector, stack them vertically and transpose into a matrix 
        # Dim of matrix = Time by Nresidues.  Residues can be inositol or protein residues.
        
        print "Saving data from", files, "into", h5file.root, table_name

        data_all = numpy.transpose(numpy.vstack(column_stack))
        atom = tables.Atom.from_dtype(data_all.dtype)
        filters = tables.Filters(complib='zlib', complevel = 5)
        h5_carray = h5file.createCArray(h5file.root, table_name, atom, data_all.shape, filters=filters)
        h5_carray[:] = data_all
        h5_carray.flush()

    def get_data_as_matrix():
        # Return the numpy representation of the matrix in memory
        return numpy.transpose(numpy.vstack(column_stack))

# Defines an analysis base class which represents a specific analysis for a system
class Analysis:
    def __init__(self, analysis_id, isomer, ratio, num_systems):
        self.analysis_id = analysis_id
        self.isomer = isomer
        self.ratio = ratio
        self.num_systems = num_systems

    def get_file_names(self, system_id):
        pass

class NonpolarAnalysisLigand:
    def get_file_names():
        pass

class NonpolarAnalysisResidue:
    def get_file_names():
        pass

class HBondAnalysisLigand:
    def __init__(self, analysis_id, isomer, ratio, num_systems, num_ligands):
        super(HBondAnalysisLigand, self).__init__(analysis_id, isomer, ratio, num_systems)
        self.num_ligands = num_ligands

    def get_file_names(self, system_id):
        data_path = os.path.join(self.analysis_name, str(system_id))
        files = [ os.path.join(data_path, str(i) + '.xvg') for i in range(1, self.num_ligands + 1) ]
        return files

class HBondAnalysisResidue:
    def __init__(self, analysis_id, isomer, ratio, num_systems, num_residues):
        super(HBondAnalysisResidue, self).__init__(analysis_id, isomer, ratio, num_systems)
        self.num_residues = num_residues

    def get_file_names(self, system_id):    
        data_path = os.path.join(self.analysis_name, str(system_id))

        # Each file here represents either a single inositol molecule or residue and the number of hbonds it froms with the aggregate
        # Form a matrix with column ordering of 1 ... 64 and store this matrix for which will be used for further post-processing
        files = [ os.path.join(data_path, str(i) + '.xvg') for i in range(0, self.num_residues) ]
        return files

class Datastore:
    def __init__(self, sys_id, isomer, ratio, analysis_name, h5file):
        self.h5file = h5file
        self.isomer = isomer
        self.ratio = ratio
        self.analysis_name = analysis_name
        self.system_id = sys_id

    def _get_table_name(self, sys_idx):
        return self.isomer + '_' + str(self.ratio) + '_' + self.analysis_name + '_' + str(self.system_id)

    def _get_tar_file_name(self):
        return "analysis_{0}_{1}_{2}.tgz".format(self.isomer, self.ratio, self.analysis_name)

    def _get_data_column_with_time(self, data):
        return [data[:, 0], data[:, 1]]

    def _get_data_column(self, data)
        return data[:, 1]

    def store_columnar_files_in_tar_to_h5(tar_obj, file_names, analysis_name, skip_header=0, analysis_name=""):
        first_file = True
        data = DataMatrix(self.h5file)
        for f in file_names:
            try:
                # TODO: Figure out the number of lines to skip
                file_data = numpy.genfromtxt(tar_obj.extractfile(f), skip_header=skip_header)
            except KeyError, e:
                print "Couldn't find", f, "in archive. Stopping read."
                return

            print file_data.shape

            if first_file:
                data.append(self._get_data_column_with_time(file_data))
            else:
                data.append(self._get_data_column(file_data))

        data.save(self._get_table_name())


class HBondDatastore(Datastore):
    def __init__(self, sys_id, isomer, ratio, analysis, h5file):
        super(HBondDatastore, self).__init__( sys_id, isomer, ratio, analysis.name(), h5file)

    def store_hbond_data(self):
        print "Munging files for", analysis.name()

        tar_file_name = get_tar_file_name(self.isomer, self.ratio, analysis.name())
        tar_file_path = os.path.join(config.data_source, tar_file_name)

        if not os.path.exists(tar_file_path):
            print "copying", tar_file_name, "to", config.data_source
            os.system("cp {0} {1}".format(tar_file_name, config.data_source))

        # TODO: Check if requires exception handling
        tar_obj = tarfile.open(tar_file_path)

        # TODO: Generalized into a function which combines a bunch of single column file into a DataMatrix object for storage into a h5file.
        for sys_idx in range(10):
            self.combine_columns_in_tar_to_data_matrix(tar_obj, analysis.get_file_names(sys_idx), skip_header=100)
            


# TODO: Refactor this into classes
def process_nonpolar_residue(h5file, isomer, ratio, analysis_type):
    print "process_nonpolar: munging files for", analysis_type


    tar_file_name = "analysis_{0}_{1}_nonpolar.tgz".format(isomer, ratio)
    tar_file_path = os.path.join(config.data_source, tar_file_name)

    if not os.path.exists(tar_file_path):
        print "copying", tar_file_name, "to", config.data_source
        os.system("cp {0} {1}".format(tar_file_name, config.data_source))
        
    # read in each of the chains and sum up the matrices
    tar = tarfile.open(tar_file_path)
    # Sum data files from all chains
    # Each file contains the inositol nonpolar contacts to a single chain in the protofibril
    # ie. for each inositol molecule is it bound to chain X, where X=1,2,3,4,5?
    # Note that the sum of all of these matrices is the total number of 
    # nonpolar contacts made by inositol to the protofibril
    for sys_idx in range(10):
        data_all = []
        read_fail = False
        for ch in range(5):
            data_file = os.path.join('nonpolar', '%(sys_idx)s_chain%(ch)s_%(analysis_type)s_np_contact' % vars() + '.dat')

            # extract the file from the tar as a file object
            print "Extracting", data_file

            try:
                member = tar.extractfile(data_file)
            except KeyError, e:
                print "Couldn't find", data_file, "in archive", "skipping ..."
                read_fail = True
                break 

            a_chain_data = numpy.genfromtxt(member, comments="#")
            if ch == 0:
                data_all.append(a_chain_data)
            else:
                data_all.append(a_chain_data[:, 1:])

        if read_fail == False:
            table_name = isomer + '_' + str(ratio) + '_' + 'nonpolar' + '_' + analysis_type + '_' + str(sys_idx)

            print "Saving data into", h5file.root, table_name

            # Saves into h5 file using pytables as a CArray object. Code straight out of documentation
            # http://readthedocs.org/docs/pytables/en/latest/usersguide/libref.html?highlight=carray#tables.CArray
            # CArray chosen for homogeneous data
            data_matrix = numpy.hstack(data_all)
            print data_matrix.shape
            
            atom = tables.Atom.from_dtype(data_matrix.dtype)
            filters = tables.Filters(complib='zlib', complevel = 5)
            h5_carray = h5file.createCArray(h5file.root, table_name, atom, data_matrix.shape, filters=filters)
            h5_carray[:] = data_matrix
            h5_carray.flush()

# TODO: Refactor this into classes        
def process_nonpolar(h5file, isomer, ratio, analysis_type):
    print "process_nonpolar: munging files for", analysis_type


    tar_file_name = "analysis_{0}_{1}_nonpolar.tgz".format(isomer, ratio)
    tar_file_path = os.path.join(config.data_source, tar_file_name)

    if not os.path.exists(tar_file_path):
        print "copying", tar_file_name, "to", config.data_source
        os.system("cp {0} {1}".format(tar_file_name, config.data_source))

    # read in each of the chains and sum up the matrices
    tar = tarfile.open(tar_file_path)
    # Sum data files from all chains
    # Each file contains the inositol nonpolar contacts to a single chain in the protofibril
    # ie. for each inositol molecule is it bound to chain X, where X=1,2,3,4,5?
    # Note that the sum of all of these matrices is the total number of 
    # nonpolar contacts made by inositol to the protofibril
    for sys_idx in range(10):
        data_all = None
        read_fail = False
        for ch in range(5):
            data_file = os.path.join('nonpolar', '%(sys_idx)s_chain%(ch)s_%(analysis_type)s_np_contact' % vars() + '.dat')

            # extract the file from the tar as a file object
            print "Extracting", data_file

            try:
                member = tar.extractfile(data_file)
            except KeyError, e:
                print "Couldn't find", data_file, "in archive", "skipping ..."
                read_fail = True
                break

            a_chain_data = numpy.genfromtxt(member, comments="#")

            if ch == 0:
                data_all = a_chain_data
            else:
                data_all[:, 1:] += a_chain_data[:, 1:]

        if read_fail == False:
            table_name = isomer + '_' + str(ratio) + '_' + 'nonpolar' + '_' + analysis_type + '_' + str(sys_idx)

            print "Saving data into", h5file.root, table_name

            # Saves into h5 file using pytables as a CArray object. Code straight out of documentation
            # http://readthedocs.org/docs/pytables/en/latest/usersguide/libref.html?highlight=carray#tables.CArray
            # CArray chosen for homogeneous data
            atom = tables.Atom.from_dtype(data_all.dtype)
            filters = tables.Filters(complib='zlib', complevel = 5)
            h5_carray = h5file.createCArray(h5file.root, table_name, atom, data_all.shape, filters=filters)
            h5_carray[:] = data_all
            h5_carray.flush()
                
                
                     
def main():
    if len(sys.argv) < 2:
        print "usage: clean.py <isomer> <ratio> <analysis>"
        sys.exit(0)

    isomer = sys.argv[1]
    ratio = int(sys.argv[2])
    analysis = sys.argv[3] 
    
    h5file = tables.openFile(isomer + "_" + str(ratio) + "_" + analysis +".h5", 'a')
    
    # reading hbonds for all inositol molecules
    if analysis == "hbonds_inositol":
        process_hbonds(h5file, isomer, ratio, "inositol")
    
    if analysis == "hbonds_residue":
        process_hbonds(h5file, isomer, ratio, "residue")
    
    if analysis == "nonpolar_inositol":
        process_nonpolar(h5file, isomer, ratio, "inositol")
    
    if analysis == "nonpolar_residue":
        process_nonpolar_residue(h5file, isomer, ratio, "residue")
    
    h5file.close()

if __name__ == '__main__':
    main()
