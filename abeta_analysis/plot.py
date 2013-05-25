#!/usr/bin/env python
import pylab
import numpy
import csv
import re
import glob
import tables
import sys

class ContactHistogram(object):
    def __init__(self, data, labels=None):
        self.data = data
        self.labels = labels

class ContactMatrix(object):
    pass

class NonpolarContactMatrix(ContactMatrix):     
    # This class represents a nonpolar contact matrix object.     
    # Data for each residue is represented by a dictionary which      
    # maps the residue name (3-letter code + residue index) to a numpy     
    # array of number which is the contacts     
    def __init__(self, h5file, table_name):
        self.h5file = h5file
        self.table_name = table_name
        self.data_table = self.h5file.getNode('/' + self.table_name).read().view(dtype=numpy.float64)
        self.contact_matrix_as_dict = {}
        self.contact_histogram = None
        self.contact_matrix = None
        self.matrix_header = "ACE0 ACE108 ACE27 ACE54 ACE81 ALA107 ALA113 ALA122 ALA134 ALA14 ALA26 ALA32 ALA41 ALA5 ALA53 ALA59 ALA68 ALA80 ALA86 ALA95 ASN11 ASN119 ASN38 ASN65 ASN92 ASP115 ASP34 ASP61 ASP7 ASP88 GLU114 GLU33 GLU6 GLU60 GLU87 ILE106 ILE123 ILE124 ILE133 ILE15 ILE16 ILE25 ILE42 ILE43 ILE52 ILE69 ILE70 ILE79 ILE96 ILE97 LEU1 LEU109 LEU126 LEU18 LEU28 LEU45 LEU55 LEU72 LEU82 LEU99 LYSH12 LYSH120 LYSH39 LYSH66 LYSH93 MET100 MET127 MET19 MET46 MET73 PHE111 PHE112 PHE3 PHE30 PHE31 PHE4 PHE57 PHE58 PHE84 PHE85 SER10 SER118 SER37 SER64 SER91 VAL101 VAL104 VAL105 VAL110 VAL116 VAL128 VAL131 VAL132 VAL2 VAL20 VAL23 VAL24 VAL29 VAL35 VAL47 VAL50 VAL51 VAL56 VAL62 VAL74 VAL77 VAL78 VAL8 VAL83 VAL89".split()

        self._prepare()

    def _prepare(self):
        matrix_dim = self.data_table.shape
        ncols = matrix_dim[1] - 1 # subtract 1 for the time column
        self.data_table_sum = numpy.average(self.data_table, axis=0)
        for col in range(ncols):
            self.contact_matrix_as_dict[ self.matrix_header[col] ] = self.data_table_sum[ col + 1 ]

        # print self.contact_matrix_as_dict

    def _compact_matrix(self, matrix):
        # Compacts the matrix. i.e. folds the matrix such that its dimensions are (rows=num_chains, column=num_residues)
        row_size = len(self.data_table[0]) / 5
        return matrix.view(dtype=numpy.float64).reshape(-1, row_size)


    def _sort_columns(self):
        data_in_order = []
        for key in sorted(self.matrix_header, key=lambda k: int(re.search('\d+', k).group(0))):
            #return data in the order of residue number and normalize by the number of atoms 
            data_column = self.contact_matrix_as_dict[key]
            data_in_order.append(data_column)

        # returns a matrix of dimension len(data_in_order) by data_column.size
        return numpy.array(data_in_order)

    def compute_contact_matrix(self):
        # Returns the contact matrix where each row corresponds to a chain
        # and each column corresponds to a residue
        # each cell is the number of nonpolar contacts to that matrix
        if self.contact_matrix is None:
            data_matrix_sorted = self._sort_columns()
            self.contact_matrix = self._compact_matrix(data_matrix_sorted)
        return self.contact_matrix
        

    def compute_contact_histogram(self):
        if self.histogram is None:
            histogram_data = numpy.average(self.compute_contact_matrix, axis=1)
            self.histogram = ContactHistogram(histogram_data, self.contact_matrix_as_dict.keys())
        return self.histogram


class HBondContactMatrix(ContactMatrix):
    def __init__(self, h5file, table_name):
        self.h5file = h5file
        self.table_name = table_name
        self.data_table = self.h5file.getNode('/' + self.table_name).read().view(dtype=numpy.float64)
        self.data_table_sum = None
        self.contact_histogram = None
        self.contact_matrix = None

    def compute_contact_matrix(self):
        self.data_table_sum = numpy.average(self.data_table[:,1:], axis=0)
        print self.data_table_sum
        row_size = self.data_table_sum.size / 5
        return self.data_table_sum.view(dtype=numpy.float64).reshape(-1, row_size)


def compute_nonpolar_matrices():
    nonpolar_analysis = "nonpolar_contacts"
    ratio_list = [15, 64]
    isomer_list = ["scyllo", "chiro", "glycerol"]

    for ratio in ratio_list:
        for isomer in isomer_list:
            h5file = tables.openFile(nonpolar_analysis + "_" + str(isomer) + "_" + str(ratio) + ".h5", 'a')
            m_total = numpy.zeros((5, 22))

            for i in range(10):
                m = NonpolarContactMatrix(h5file, "%(isomer)s_%(ratio)s_residue_nonpolar_contacts_%(i)d" % vars())
                matrix = m.compute_contact_matrix()
                m_total += m.compute_contact_matrix()

            m_total = m_total / 10.0
            print "calculating for", isomer, ratio
            numpy.savetxt("%(isomer)s_%(ratio)s_contact_matrix.txt" % vars(), m_total, fmt="%.2f", delimiter=' ')
            histogram = numpy.average(m_total, axis=0)
            numpy.savetxt("%(isomer)s_%(ratio)s_histogram.txt" % vars(), histogram, fmt="%.2f", delimiter=' ')


def compute_hbond_matrices():
    analysis = "hbonds"
    ratio_list = [15, 64]
    isomer_list = ["scyllo", "chiro", "glycerol"]

    for ratio in ratio_list:
        for isomer in isomer_list:
            h5file_name = analysis + "_" + str(isomer) + "_" + str(ratio) + ".h5"
            print h5file_name
            h5file = tables.openFile(h5file_name, 'a')
            m_total = numpy.zeros((5, 26))

            for i in range(10):
                m = HBondContactMatrix(h5file, "%(isomer)s_%(ratio)s_residue_hbonds_%(i)d" % vars())
                matrix = m.compute_contact_matrix()
                m_total += m.compute_contact_matrix()

            m_total = m_total / 10.0
            print "calculating for", isomer, ratio
            numpy.savetxt("%(isomer)s_%(ratio)s_contact_matrix.txt" % vars(), m_total, fmt="%.2f", delimiter=' ')
            histogram = numpy.average(m_total, axis=0)
            numpy.savetxt("%(isomer)s_%(ratio)s_histogram.txt" % vars(), histogram, fmt="%.2f", delimiter=' ')


if __name__ == '__main__':
    compute_hbond_matrices()
