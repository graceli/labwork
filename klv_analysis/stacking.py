import tables
import numpy

import plot_and_save2hdf5 as myh5


# Returns the number of stacked and bound phes for a given data vector (corresponds to a result from a trajectory frame)
def match_phe_binding(residue_dict, stacking_dict):
    for keys in residue_dict:
        if residue_dict[key] > 0 and stacking_dict[key] > 0:
            stacked += 1
        
        if residue_dict[key] > 0:
            bound += 1
    
    return stacked,bound
    
def beta_stacking(h5file, type, system_indices = []):
    phe_header_residue = "PHE103 PHE104 PHE112 PHE113 PHE121 PHE122 PHE13 PHE130 PHE131 PHE139 PHE14 PHE140 PHE22 PHE23 PHE31 PHE32 PHE4 PHE40 PHE41 PHE49 PHE5 PHE50 PHE58 PHE59 PHE67 PHE68 PHE76 PHE77 PHE85 PHE86 PHE94 PHE95".split()
    phe_header_stacking = "PHE4 PHE5 PHE13 PHE14 PHE22 PHE23 PHE31 PHE32 PHE40 PHE41 PHE49 PHE50 PHE58 PHE59 PHE67 PHE68 PHE76 PHE77 PHE85 PHE86 PHE94 PHE95 PHE103 PHE104 PHE112 PHE113 PHE121 PHE122 PHE130 PHE131 PHE139 PHE140".split()
    for i in system_indices:
        stacked_system_total = 0
        bound_system_total = 0

        residue_file = '/nonpolar_residue_dt1/scyllo_t%(i)d_per_residue_contacts.dat' % vars()
        phe_stacking_file = '/stacking/scyllo_sys%(i)d_per_phe_stacking.dat' % vars()
        
        residue_matrix = myh5.getTableAsMatrix(h5file, residue_file, dtype=numpy.float64)
        phe_stacking  = myh5.getTableAsMatrix(h5file, phe_stacking_file, dtype=numpy.float64)
        
        nrows, ncols = residue_matrix.shape
        for row in range(0, nrows):
            residue_data_list = numpy.array(residue_matrix[row, 64:96].T).flatten().tolist()
            stacking_data_list = numpy.array(phe_stacking[row, :].T).flatten().tolist()
            
            residue_data_dict = dict(zip(phe_header_residue, residue_data_list))
            stacking_data_dict = dict(zip(phe_stacking, stacking_data_list))
            
            stacked, bound = match_phe_binding(residue_data_dict, stacking_data_dict)
            
            stacked_system_total += stacked
            bound_system_total += bound

        print stacked_system_total, bound_system_total, stacked_system_total / float(bound_system_total)


def oligomer_stacking(h5file, type, system_indices = []):
    phe_header_residue = "PHE13 PHE14 PHE22 PHE23 PHE31 PHE32 PHE4 PHE5".split()
    phe_header_stacking = "PHE4 PHE5 PHE13 PHE14 PHE22 PHE23 PHE31 PHE32".split()
    for i in system_indices:
        stacked_system_total = 0
        bound_system_total = 0
        residue_file = '/nonpolar_residue/scyllo_sys%(i)d_%(type)d_whole_nosol_0-200ns_per_residue_contact.dat' % vars()

        # /stacking/scyllo_sys9_per_phe_stacking.dat
        phe_stacking_file = '/stacking/scyllo_sys%(i)d_per_phe_stacking.dat' % vars()

        residue_matrix = myh5.getTableAsMatrix(h5file, residue_file, dtype=numpy.float64)
        phe_stacking  = myh5.getTableAsMatrix(h5file, phe_stacking_file, dtype=numpy.float64)

        nrows,ncols = residue_matrix.shape
        for row in range(0, nrows):
            # Grab a row from a numpy matrix and converts it to a list
            # Note that numpy matrix slicing range is not inclusive at the higher index
            # 19:24 gets the columns [19,23]
            residue_data_list = numpy.array(residue_matrix[row, 16:24].T).flatten().tolist()
            stacking_data_list = numpy.array(phe_stacking[row, :].T).flatten().tolist()
            
            residue_data_dict = dict(zip(phe_header_residue, residue_data_list))
            stacking_data_dict = dict(zip(phe_stacking, stacking_data_list))
            
            stacked,bound = match_phe_binding(residue_data_dict, stacking_data_dict)

            stacked_system_total += stacked
            bound_system_total += bound

        print stacked_system_total, bound_system_total, stacked_system_total / float(bound_system_total)


# TODO: refactor this to work for both ratios for the monomer system
def monomer_stacking_2to1(h5file):
    # h5file = tables.openFile('monomer_2to1.h5')
        
    for i in range(0, 6):
        residue_file = '/nonpolar_residue/scyllo_sys%(i)d_mon_2to1_per_residue_contacts.dat' % vars()

        # /stacking/scyllo_sys9_per_phe_stacking.dat
        phe_stacking_file = '/stacking/scyllo_sys%(i)d_per_phe_stacking.dat' % vars()

        residue_matrix = myh5.getTableAsMatrix(h5file, residue_file, dtype=numpy.float64)
        phe_stacking  = myh5.getTableAsMatrix(h5file, phe_stacking_file, dtype=numpy.float64)

        nrows, ncols = residue_matrix.shape
        print residue_matrix.shape, phe_stacking.shape
        bound = 0.0
        stacked = 0.0
        for k in range(0, nrows):
            if residue_matrix[k][5] > 0: 
                bound = bound + 1

            if residue_matrix[k][6] > 0: 
                bound = bound + 1

            if phe_stacking[k][1] > 0:
                stacked = stacked + 1

            if phe_stacking[k][2] > 0:
                stacked = stacked + 1

        print stacked, bound, stacked / float(bound)