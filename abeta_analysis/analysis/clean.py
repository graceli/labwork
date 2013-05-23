import os 
import tarfile
import logging
import re

import numpy 
import tables

import config 

# Note: all my operations involves the following
# can I build a system (or is there an existing one) to suit my needs?
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
        
        print "Saving data into", self.h5file.root, table_name

        atom = tables.Atom.from_dtype(data_matrix.dtype)
        filters = tables.Filters(complib='zlib', complevel=5)
        h5_carray = self.h5file.createCArray(self.h5file.root, table_name, atom, data_matrix.shape, filters=filters)
        h5_carray[:] = data_matrix
        h5_carray.flush()


class DataMatrixFromColumnsBuilder:
    def __init__(self):
        self.column_stack = []

    def clear(self):
        self.column_stack = []

    def append(self, array):
        # Appends a column to the data matrix
        self.column_stack.append(array)

    def get_data_as_numpy_matrix(self):
        # Return the numpy representation of the matrix in memory
        return numpy.transpose(numpy.vstack(self.column_stack))


# Defines an analysis base class which represents a specific analysis for a system
class Analysis(object):
    def __init__(self, name, isomer, ratio, num_systems, num_residues):
        self.name = name
        self.isomer = isomer
        self.ratio = ratio
        self.num_systems = num_systems
        self.num_residues = num_residues

    def _get_table_name(self, file_name):
        m = re.search("([a-zA-Z_]*)(\d+)([a-zA-Z_]*)(\d+)", file_name)
        # Get the matching string at the 3rd position
        system_id = int(m.group(2))

        return self.isomer + '_' + str(self.ratio) + '_' + self.name + '_' + str(system_id)

    def _get_tar_file_name(self):
        return "analysis_{0}_{1}_{2}.tgz".format(self.name, self.isomer, self.ratio)

    def get_file_names(self, system_id):
        pass


class HBondAnalysisLigand(Analysis):
    def __init__(self, name, isomer, ratio, num_systems, num_ligands):
        super(HBondAnalysisLigand, self).__init__(name, isomer, ratio, num_systems, num_ligands)

    def _get_table_name(self, file_name):
        m = re.search("([a-zA-Z_]*)(\d+)([a-zA-Z_]*)(\d+)", file_name)
        # Get the matching string at the 3rd position
        system_id = int(m.group(4))

        return self.isomer + '_' + str(self.ratio) + '_inositol_' + self.name + '_' + str(system_id)

    def _get_tar_file_name(self):
        return "analysis_{0}_{1}_{2}_by_ligand.tgz".format(self.name, self.isomer, self.ratio)

    def get_file_names(self):
        files = [ 'ab_' + self.isomer + '_' + str(self.ratio) + '_' + str(idx) + '_ins' + str(lig) + '.xvg' 
        for idx in range(1, self.num_systems) for lig in range(1, self.num_residues + 1) ]

        print files
        
        return files


class HBondAnalysisResidue(Analysis):
    def __init__(self, name, isomer, ratio, num_systems, num_residues):
        super(HBondAnalysisResidue, self).__init__(name, isomer, ratio, num_systems, num_residues)

    def _get_table_name(self, file_name):
        m = re.search("([a-zA-Z_]*)(\d+)([a-zA-Z_]*)(\d+)", file_name)
        # Get the matching string at the 3rd position
        system_id = int(m.group(4))

        return self.isomer + '_' + str(self.ratio) + '_residue_' + self.name + '_' + str(system_id)

    def _get_tar_file_name(self):
        return "analysis_{0}_{1}_{2}_by_residue.tgz".format(self.name, self.isomer, self.ratio)

    def get_file_names(self):    
        files = [ 'ab_' + self.isomer + '_' + str(self.ratio) + '_' + str(idx) + '_residue' + str(lig) + '.xvg' 
        for idx in range(1, self.num_systems+1) for lig in range(0, self.num_residues) ]
        print files
        return files


class NonpolarAnalysisLigand(Analysis):
    def __init__(self, name, isomer, ratio, num_systems, num_ligands):
        super(NonpolarAnalysisLigand, self).__init__(name, isomer, ratio, num_systems, num_ligands)

    def _get_table_name(self, file_name):
        m = re.search("([a-zA-Z_]*)(\d+)([a-zA-Z_]*)(\d+)", file_name)
        # Get the matching string at the 3rd position
        system_id = int(m.group(4))

        return self.isomer + '_' + str(self.ratio) + '_inositol_' + self.name + '_' + str(system_id)

    def get_file_names(self):
        file_prefix = 'ab_' + self.isomer + '_' + str(self.ratio) + '_'

        # Each file here represents either a single inositol molecule or residue and the number of hbonds it froms with the aggregate
        # Form a matrix with column ordering of 1 ... 64 and store this matrix for which will be used for further post-processing
        files = [ file_prefix + str(i) + '_inositol_np_contact' + '.dat' for i in range(0, self.num_systems + 1) ]
        print files

        return files


class NonpolarAnalysisResidue(Analysis):
    def __init__(self, name, isomer, ratio, num_systems, num_residues):
        super(NonpolarAnalysisResidue, self).__init__(name, isomer, ratio, num_systems, num_residues)

    def _get_table_name(self, file_name):
        m = re.search("([a-zA-Z_]*)(\d+)([a-zA-Z_]*)(\d+)", file_name)
        # Get the matching string at the 3rd position
        system_id = int(m.group(4))

        return self.isomer + '_' + str(self.ratio) + '_residue_' + self.name + '_' + str(system_id)

    def get_file_names(self):
        file_prefix = 'ab_' + self.isomer + '_' + str(self.ratio) + '_'
        files = [ file_prefix + str(i) + '_residue_np_contact' + '.dat' for i in range(0, self.num_systems + 1) ]
        print files

        return files


class Datastore(object):
    def __init__(self, analysis, h5file):
        self.h5file = h5file
        self.analysis = analysis

    def store_columnar_files_in_tar_to_h5(self, tar_obj, file_names, skip_header=0):
        try:
            first_file = True
            data_builder = DataMatrixFromColumnsBuilder()
            count = 0
            for f in file_names:
                print "Storing file", f
                # TODO: Figure out the number of lines to skip
                file_data = numpy.genfromtxt(tar_obj.extractfile(f), skip_header=skip_header)
        
                print file_data.shape

                if first_file:
                    data_builder.append(self._get_data_column_with_time(file_data))
                    first_file = False
                else:
                    data_builder.append(self._get_data_column(file_data))

                count += 1

                # Need to save when all of the residues (i.e. Nresidue number of files) are added to the
                # matrix
                if count == self.analysis.num_residues:
                    data_matrix = DataMatrix(self.h5file)
                    data_matrix.save(data_builder.get_data_as_numpy_matrix(), self.analysis._get_table_name(f))
                    # TODO: Basically the count is a way to detect whether we've changed system indices. We have a System object where it is a collection of N systems with the same properties. I think the most appropriate / error free way to do this is to refactor it into objects which is iterated through by the number of systems, and returns all the file names mapping to that system
                    count = 0
                    first_file = True
                    data_builder.clear()
        except KeyError:
            # TODO: Should refactor out these error messages
            print "Couldn't find", f, "in archive. Nothing was stored. Stopping read."
            return

    def store_file_as_data_matrix(self, tar_obj, file_names, skip_header=0):
        for f in file_names:
            logging.debug("Storing file %s", f)
            try:
                file_data = numpy.genfromtxt(tar_obj.extractfile(f), skip_header=skip_header)
                data = DataMatrix(self.h5file)
                table_name = self.analysis._get_table_name(f)
                data.save(file_data, table_name)
                logging.debug("Storing file %s into table %s", f, table_name)
            except KeyError:
                print "Couldn't find", f, "in archive. Nothing was stored. Skipping"
                # NOTE: Files containing data matrices are self-contained because they constitute their own table

    # TODO: Re-implement table naming convention
    # def _get_table_name(self, analysis_file_name):
    #     name, ext = os.path.splitext(analysis_file_name)
    #     return name

    def _get_tar_obj(self):
        tar_file_name = self.analysis._get_tar_file_name()
        tar_file_path = os.path.join(os.environ['PWD'] + '/' + self.analysis.name, tar_file_name)

        if not os.path.exists(tar_file_path):
            # logging.debug("copying %s to %s", tar_file_name, config.data_source)
            # subprocess.check_call("cp {0} {1}".format(tar_file_name, config.data_source), shell=True)
            logging.debug("%s was not found", tar_file_path)

        try:
            logging.debug("Open tarfile %s", tar_file_path)
            return tarfile.open(tar_file_path)
        except tarfile.CompressionError:
            logging.error("Problem reading tar file %s", tar_file_name)
        except tarfile.ReadError:
            logging.error("Problem reading tar file %s", tar_file_name)
        
        return None

    def _get_data_column_with_time(self, data):
        return [data[:, 0], data[:, 1]]

    def _get_data_column(self, data):
        return data[:, 1]


class HBondDatastore(Datastore):
    def __init__(self, analysis, h5file):
        super(HBondDatastore, self).__init__(analysis, h5file)

    def store_hbond_data(self, skip_header=0):
        print "Munging files for", self.analysis.name
        tar_obj = self._get_tar_obj()
        self.store_columnar_files_in_tar_to_h5(tar_obj, self.analysis.get_file_names(), skip_header=skip_header)


class NonpolarDataStore(Datastore):
    def __init__(self, analysis, h5file):
        super(NonpolarDataStore, self).__init__(analysis, h5file)

    def store_nonpolar_data(self, skip_header=0):
        print "Munging files for", self.analysis.name
        tar_obj = self._get_tar_obj()
        if tar_obj:
            self.store_file_as_data_matrix(tar_obj, self.analysis.get_file_names(), skip_header=skip_header)
        else:
            print "tar file could not be opened"


def main():
    ratio = 64
    isomer = "scyllo"
    num_systems = 10
    num_ligands = 64
    num_residues = 130
    analysis_name = "hbonds"
    nonpolar_analysis_name = "nonpolar_contacts"
    # analysis_nonpolar_contacts_scyllo_64.tgz
    h5file = tables.openFile(analysis_name + "_" + str(isomer) + "_" + str(ratio) + ".h5", 'a')

    # hb_ligand = HBondAnalysisLigand(analysis_name, isomer, ratio, num_systems, num_ligands)
    # hb_ligand_store = HBondDatastore(hb_ligand, h5file)
    # hb_ligand_store.store_hbond_data(skip_header=20)

    # hb_residue = HBondAnalysisResidue(analysis_name, isomer, ratio, num_systems, num_residues)
    # hb_residue_store = HBondDatastore(hb_residue, h5file)
    # hb_residue_store.store_hbond_data(skip_header=20)

    np_ligand = NonpolarAnalysisLigand(nonpolar_analysis_name, isomer, ratio, num_systems, num_ligands)
    np_ligand_store = NonpolarDataStore(np_ligand, h5file)
    np_ligand_store.store_nonpolar_data()

    np_residue = NonpolarAnalysisResidue(nonpolar_analysis_name, isomer, ratio, num_systems, num_residues)
    np_residue_store = NonpolarDataStore(np_residue, h5file)
    np_residue_store.store_nonpolar_data()

if __name__ == '__main__':
    main()

