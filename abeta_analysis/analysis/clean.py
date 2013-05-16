import subprocess
import os
import tarfile
import sys
import numpy
import tables
import config
import logging


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

logging.basicConfig(filename='analysis.log', format='%(levelname)s:%(message)s', level=logging.DEBUG)

class DataMatrix:
    def __init__(self, h5file):
        self.h5file = h5file

    def save(self, data, table_name):
        self._save_to_h5(data, table_name)

    def _save_to_h5(self, data_matrix, table_name):
        # Each item in list is a row vector, stack them vertically and transpose into a matrix 
        # Dim of matrix = Time by Nresidues.  Residues can be inositol or protein residues.
        
        print "Saving data from", files, "into", h5file.root, table_name

        atom = tables.Atom.from_dtype(data_matrix.dtype)
        filters = tables.Filters(complib='zlib', complevel=5)
        h5_carray = h5file.createCArray(h5file.root, table_name, atom, data_matrix.shape, filters=filters)
        h5_carray[:] = data_matrix
        h5_carray.flush()


class DataMatrixFromColumnsBuilder:
    def __init__(self):
        self.column_stack = []

    def append(self, array):
        # Appends a column to the data matrix
        self.column_stack.append(array)

    def get_data_as_numpy_matrix(self):
        # Return the numpy representation of the matrix in memory
        if self.data is not None:
            return self.data
        return numpy.transpose(numpy.vstack(self.column_stack))


# Defines an analysis base class which represents a specific analysis for a system
class Analysis:
    def __init__(self, analysis_id, isomer, ratio, num_systems):
        self.analysis_id = analysis_id
        self.isomer = isomer
        self.ratio = ratio
        self.num_systems = num_systems

    def get_file_names(self, system_id):
        pass


class HBondAnalysisLigand(Analysis):
    def __init__(self, analysis_id, isomer, ratio, num_systems, num_ligands):
        super(HBondAnalysisLigand, self).__init__(analysis_id, isomer, ratio, num_systems)
        self.num_ligands = num_ligands

    def get_file_names(self, system_id):
        data_path = os.path.join(self.analysis_name, str(system_id))
        files = [ os.path.join(data_path, str(i) + '.xvg') for i in range(1, self.num_ligands + 1) ]
        return files


class HBondAnalysisResidue(Analysis):
    def __init__(self, analysis_id, isomer, ratio, num_systems, num_residues):
        super(HBondAnalysisResidue, self).__init__(analysis_id, isomer, ratio, num_systems)
        self.num_residues = num_residues

    def get_file_names(self, system_id):    
        data_path = os.path.join(self.analysis_name, str(system_id))

        # Each file here represents either a single inositol molecule or residue and the number of hbonds it froms with the aggregate
        # Form a matrix with column ordering of 1 ... 64 and store this matrix for which will be used for further post-processing
        files = [ os.path.join(data_path, str(i) + '.xvg') for i in range(0, self.num_residues) ]
        return files


class NonpolarAnalysisLigand(Analysis):
    def __init__(self, analysis_id, isomer, ratio, num_systems, num_ligands):
        super(HBondAnalysisLigand, self).__init__(analysis_id, isomer, ratio, num_systems)
        self.num_ligands = num_ligands

    def get_file_names(self):
        data_path = os.path.join(self.analysis_name, str(system_id))
        files = [ os.path.join(data_path, str(i) + '.xvg') for i in range(1, self.num_ligands + 1) ]
        return files


class NonpolarAnalysisResidue(Analysis):
    def __init__(self, analysis_id, isomer, ratio, num_systems, num_residues):
        super(HBondAnalysisResidue, self).__init__(analysis_id, isomer, ratio, num_systems)
        self.num_residues = num_residues

    def get_file_names(self):
        data_path = os.path.join(self.analysis_name, str(system_id))

        # Each file here represents either a single inositol molecule or residue and the number of hbonds it froms with the aggregate
        # Form a matrix with column ordering of 1 ... 64 and store this matrix for which will be used for further post-processing
        files = [ os.path.join(data_path, str(i) + '.xvg') for i in range(0, self.num_residues) ]
        return files


class Datastore:
    def __init__(self, analysis, h5file):
        self.h5file = h5file
        self.analysis = analysis

    def store_columnar_files_in_tar_to_h5(self, tar_obj, file_names, skip_header=0):
        first_file = True
        data_builder = DataMatrixFromColumnsBuilder(self.h5file)
        for f in file_names:
            try:
                # TODO: Figure out the number of lines to skip
                file_data = numpy.genfromtxt(tar_obj.extractfile(f), skip_header=skip_header)
            except KeyError as e:
                print "Couldn't find", f, "in archive. Stopping read."
                return

            print file_data.shape

            if first_file:
                data_builder.append(self._get_data_column_with_time(file_data))
            else:
                data_builder.append(self._get_data_column(file_data))

        data_matrix = DataMatrix()
        data_matrix.save(data_builder.get_data_as_matrix(), self._get_table_name())


    def store_file_as_data_matrix(self, tar_obj, file_names, skip_header=0):
        for f in file_names:
            logging.debug("Storing file %s", f)

            try:
                file_data = numpy.genfromtxt(tar_obj.extractfile(f), skip_header=skip_header)
            except KeyError as e:
                print "Couldn't find", f, "in archive. Stopping read."
                return

            data = DataMatrix(self.h5file)
            data.save(file_data, self._get_table_name())


    def _get_table_name(self, sys_idx):
        return self.isomer + '_' + str(self.analysis.ratio) + '_' + self.analysis.name + '_' + str(self.analysis.system_id)

    def _get_tar_file_name(self):
        return "analysis_{0}_{1}_{2}.tgz".format(self.analysis.name, self.analysis.isomer, self.analysis.ratio)

    def _get_tar_obj(self):
        tar_file_name = self._get_tar_file_name()
        tar_file_path = os.path.join(config.data_source, tar_file_name)

        if not os.path.exists(tar_file_path):
            logging.debug("copying %s to %s", tar_file_name, config.data_source)
            subprocess.check_call("cp {0} {1}".format(tar_file_name, config.data_source), shell=True)
        else:
            logging.debug("%s was not found", tar_file_path)
            return

        try:
            tar_obj = tarfile.open(tar_file_path)
        except CompressError:
            logging.error("Problem reading tar file %s", tar_file_name)
            return
        except ReadError:
            logging.error("Problem reading tar file %s", tar_file_name)
            return

    def _get_data_column_with_time(self, data):
        return [data[:, 0], data[:, 1]]

    def _get_data_column(self, data):
        return data[:, 1]


class HBondDatastore(Datastore):
    def __init__(self, analysis, h5file):
        super(HBondDatastore, self).__init__(analysis, h5file)

    def store_hbond_data(self):
        print "Munging files for", self.analysis.name
        tar_obj = self._get_tar_obj()

        for sys_idx in range(num_systems):
            self.store_columnar_files_in_tar_to_h5(tar_obj, self.analysis.get_file_names(sys_idx), skip_header=100)


class NonpolarDataStore(Datastore):
    def __init__(self, analysis, h5file):
        super(NonpolarDataStore, self).__init__(analysis, h5file)

    def store_nonpolar_data(self):
        print "Munging files for", self.analysis.name
        tar_obj = self._get_tar_obj()
        for sys_idx in range(num_systems):
            self.store_file_as_data_matrix(tar_obj, self.analysis.get_file_names(sys_idx), skip_header=100)


def main():
    index = 0
    ratio = 15
    isomer = "scyllo"
    num_systems = 10
    num_ligands = 15
    analysis_name = "nonpolar_contacts"
    # analysis_nonpolar_contacts_scyllo_64.tgz
    h5file = tables.openFile(analysis_name + "_" + str(isomer) + "_" + str(ratio) + ".h5", 'a')
    for sys_id in range(num_systems):
        print sys_id
        nonpolar_analysis = NonpolarAnalysisLigand(analysis_name, isomer, ratio, num_systems, num_ligands)
        nonpolar = NonpolarDataStore(nonpolar_analysis, h5file)
        nonpolar.store_nonpolar_data()


def '__name__' == '__main__':

