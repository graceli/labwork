import tables
import numpy
import sys
import csv
import os

import plot_and_save2hdf5 as myh5


# Returns the number of stacked and bound phes for a given data vector (corresponds to a result from a trajectory frame)
def match_phe_binding(residue_dict, stacking_dict):
    stacked = 0
    bound = 0
    stacked_bound = 0
    for key in residue_dict:
        if stacking_dict[key] > 0:
            stacked += 1
       
        if stacking_dict[key] > 0 and residue_dict[key] > 0:
            stacked_bound += 1 

        if residue_dict[key] > 0:
            bound += 1
    
    return stacked,bound,stacked_bound
    
def beta_stacking(h5file, type, system_indices = [], tag="15", file_path='/stacking'):
    writer = csv.writer(open(type + "_" + tag + "_beta_stacking.csv", 'wb'))
    writer.writerow(["stacked", "bound", "stacked+bound", "stacked/bound"])

    phe_header_residue = "PHE103 PHE104 PHE112 PHE113 PHE121 PHE122 PHE13 PHE130 PHE131 PHE139 PHE14 PHE140 PHE22 PHE23 PHE31 PHE32 PHE4 PHE40 PHE41 PHE49 PHE5 PHE50 PHE58 PHE59 PHE67 PHE68 PHE76 PHE77 PHE85 PHE86 PHE94 PHE95".split()
    phe_header_stacking = "PHE4 PHE5 PHE13 PHE14 PHE22 PHE23 PHE31 PHE32 PHE40 PHE41 PHE49 PHE50 PHE58 PHE59 PHE67 PHE68 PHE76 PHE77 PHE85 PHE86 PHE94 PHE95 PHE103 PHE104 PHE112 PHE113 PHE121 PHE122 PHE130 PHE131 PHE139 PHE140".split()
    for i in system_indices:
        stacked_system_total = 0
        bound_system_total = 0
        stacked_bound_system_total = 0
        residue_file = ""
        if type == 'low':
            residue_file = '/nonpolar_residue_dt1/scyllo_t%(i)d_per_residue_contacts.dat' % vars()
        elif type == 'high':
            residue_file = '/nonpolar_residue/scyllo_t%(i)d_per_residue_contacts.dat' % vars()
        else:
            print "system type", type, "is not recognize"
            sys.exit()

        phe_stacking_file = os.path.join(file_path, 'scyllo_sys%(i)d_per_phe_stacking.dat' % vars())

        print residue_file, phe_stacking_file
 
        residue_matrix = myh5.getTableAsMatrix(h5file, residue_file, dtype=numpy.float64)
        phe_stacking  = myh5.getTableAsMatrix(h5file, phe_stacking_file, dtype=numpy.float64)
        
        print residue_matrix.shape, phe_stacking.shape
        assert residue_matrix.shape[0] == phe_stacking.shape[0], "residue matrix and phe_stacking matrix must have the same number of rows" 
  
        nrows, ncols = residue_matrix.shape
        for row in range(0, nrows):
            residue_data_list = numpy.array(residue_matrix[row, 65:97].T).flatten().tolist()
            stacking_data_list = numpy.array(phe_stacking[row, 1:].T).flatten().tolist()
            
            residue_data_dict = dict(zip(phe_header_residue, residue_data_list))
            stacking_data_dict = dict(zip(phe_header_stacking, stacking_data_list))
            
            stacked,bound,stacked_bound = match_phe_binding(residue_data_dict, stacking_data_dict)
            
            stacked_system_total += stacked
            bound_system_total += bound
            stacked_bound_system_total += stacked_bound

        writer.writerow([stacked_system_total, bound_system_total, stacked_bound_system_total, stacked_system_total / float(bound_system_total)])


def oligomer_stacking(h5file, type, system_indices = [], tag="15", file_path='/stacking'):
    phe_header_residue = "PHE13 PHE14 PHE22 PHE23 PHE31 PHE32 PHE4 PHE5".split()
    phe_header_stacking = "PHE4 PHE5 PHE13 PHE14 PHE22 PHE23 PHE31 PHE32".split()

    writer = csv.writer(open(type + "_" + tag + "_oligomer_stacking.csv", 'wb'))
    writer.writerow(["stacked", "bound", "stacked+bound", "stacked/bound"])
    for i in system_indices:
        stacked_system_total = 0
        bound_system_total = 0
        stacked_bound_system_total = 0
        residue_file = ""
        if type == "15to4":
            residue_file = '/nonpolar_residue/scyllo_sys%(i)d_%(type)s_whole_nosol_0-200ns_per_residue_contact.dat' % vars()
        elif type == "45to4":            
           residue_file = '/nonpolar_residue/scyllo_sys%(i)d_%(type)s_whole_nosol_0-200_per_residue_contact.dat' % vars()
        else:
            print "unrecognized system type", type
            sys.exit()

        # /stacking/scyllo_sys9_per_phe_stacking.dat
        phe_stacking_file = os.path.join(file_path, 'scyllo_sys%(i)d_per_phe_stacking.dat' % vars())

        residue_matrix = myh5.getTableAsMatrix(h5file, residue_file, dtype=numpy.float64)
        phe_stacking  = myh5.getTableAsMatrix(h5file, phe_stacking_file, dtype=numpy.float64)

        if residue_matrix is None or phe_stacking is None:
            print residue_file, phe_stacking_file, "does not exist"
            continue

        nrows,ncols = residue_matrix.shape
        for row in range(0, nrows):
            # Grab a row from a numpy matrix and converts it to a list
            # Note that numpy matrix slicing range is not inclusive at the higher index
            # Note that the column indices accounts for the first column as being time
            residue_data_list = numpy.array(residue_matrix[row, 17:25].T).flatten().tolist()
            stacking_data_list = numpy.array(phe_stacking[row, 1:].T).flatten().tolist()
            
            residue_data_dict = dict(zip(phe_header_residue, residue_data_list))
            stacking_data_dict = dict(zip(phe_header_stacking, stacking_data_list))
            
            stacked,bound,stacked_bound = match_phe_binding(residue_data_dict, stacking_data_dict)

            stacked_system_total += stacked
            bound_system_total += bound
            stacked_bound_system_total += stacked_bound

        writer.writerow([stacked_system_total, bound_system_total, stacked_bound_system_total, stacked_system_total / float(bound_system_total)])

def monomer_stacking(h5file, ratio, system_indices, tag="15", file_path='/stacking'):
    writer = csv.writer(open(ratio + "_" + tag + "_monomer_stacking.csv", 'wb'))
    header = ["stacked","bound", "stacked+bound", "stacked/bound"]

    writer.writerow(header)
    for i in system_indices:
        residue_file = ""
        phe_stacking_file = ""
        if ratio == "2to1":
            residue_file = '/nonpolar_residue/scyllo_sys%(i)d_mon_2to1_per_residue_contacts.dat' % vars()
            phe_stacking_file = os.path.join(file_path, 'scyllo_sys%(i)d_per_phe_stacking.dat') % vars()
            print residue_file, phe_stacking_file
        elif ratio == "15to1":
            residue_file = '/nonpolar/scyllo_sys%(i)d_per_residue_contacts.dat' % vars()
            phe_stacking_file = os.path.join(file_path, 'scyllo_sys%(i)d_per_phe_stacking.dat') % vars()
            print residue_file, phe_stacking_file
        else:
            # TODO: Throw a custom exception here
            print "ratio ", ratio, "is not recognized"
            sys.exit()
            
        residue_matrix = myh5.getTableAsMatrix(h5file, residue_file, dtype=numpy.float64)
        phe_stacking  = myh5.getTableAsMatrix(h5file, phe_stacking_file, dtype=numpy.float64)

        if residue_matrix is None or phe_stacking is None:
            print residue_file, "or", phe_stacking_file, "does not exist"
            continue

        print residue_matrix.shape, phe_stacking.shape
        assert residue_matrix.shape[0] == phe_stacking.shape[0], "Residue and phe stacking matrices must have the same number of lines"

        nrows, ncols = residue_matrix.shape

        bound = 0.0
        stacked = 0.0
        stacked_bound = 0.0
        for k in range(0, nrows):
            if residue_matrix[k][5] > 0: 
                bound = bound + 1

            if residue_matrix[k][6] > 0: 
                bound = bound + 1

            if phe_stacking[k][1] > 0:
                stacked = stacked + 1

            if phe_stacking[k][2] > 0:
                stacked = stacked + 1

            # this is for a sanity check
            if residue_matrix[k][5] > 0 and phe_stacking[k][1] > 0:
                stacked_bound += 1
            
            if residue_matrix[k][6] > 0 and phe_stacking[k][2] > 0:
                stacked_bound += 1

        writer.writerow([stacked, bound, stacked_bound, stacked / float(bound)])

