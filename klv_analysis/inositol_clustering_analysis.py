import csv
import tables
import sys
import os
import numpy

import plot_and_save2hdf5 as myh5

# Returns true if the inositol cluster is bound to protein or not
def _residue_ids_to_indices(inositol_residue_ids):
    # The offset is the starting residue id of inositol molecules in the trajectory - eg. 145 for beta in the gro file (144 is the id of the first inositol molecule in the gromacs analysis framework). The indices in the inositol contacts data file is indexed by the column number. ie. 0th column is the 0th inositol corresponding to inositol residue id of 144.
    # TODO: For now hardcode this offset to the value for the beta system
    offset = 144
    inositol_cluster_indices = [int(x) - offset for x in inositol_residue_ids]
    return inositol_cluster_indices


def _contacts_for_indices_in_cluster(inositol_contacts_array, inositol_cluster_indices):
    return inositol_contacts_array[inositol_cluster_indices]

def _cluster_bound_to_protein(inositol_in_cluster_contacts):
    num_contacts_with_cluster =  numpy.sum(inositol_in_cluster_contacts)
    return num_contacts_with_cluster > 0
    
# for iso in isomerList:
#      for sys in system_indices[iso]:

def compute_inositol_ub_b_cluster_size_histo(h5file, clust_info_path, iso, sys, tag=""):
    polar_contacts_file = os.path.join('/polar', "%(iso)s_t%(sys)d_inos_total.dat" % vars())
    nonpolar_contacts_file = os.path.join('/nonpolar_residue', "%(iso)s_t%(sys)d_per_inositol_contacts.dat" % vars())

#/nonpolar_residue_dt1/chiro_t1_per_inositol_contacts.dat

    polar_contacts = myh5.getTableAsMatrix(h5file, polar_contacts_file, dtype=numpy.float64)
    nonpolar_contacts = myh5.getTableAsMatrix(h5file, nonpolar_contacts_file, dtype=numpy.float64)

    clust_info_csv = os.path.join(clust_info_path, '%(iso)s_high_conc_t%(sys)d_all_nosol_whole_clust_info.dat' % vars())
   
    contacts_matrix = polar_contacts[:,1:] + nonpolar_contacts[:,1:]

    print contacts_matrix.shape

    # Parse the csv inositol clusters data
    bound_sizes_list = []
    unbound_sizes_list = []
    
    with open(clust_info_csv, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            time = row[0]
            clustered = row[1]
            inositol_ids = row[2].split(' ')
           
            # print time, clustered, inositol_ids
 
            inositol_indices = _residue_ids_to_indices(inositol_ids)
            try:
                inos_in_cluster_contacts = _contacts_for_indices_in_cluster(contacts_matrix[int(time) / 2], inositol_indices)
            except IndexError:
                print "index error encountered at time=", time
                break
            
            if clustered == "yes":
                bound = _cluster_bound_to_protein(inos_in_cluster_contacts)
                if bound:
                    bound_sizes_list.append(len(inositol_indices))
                else:
                    unbound_sizes_list.append(len(inositol_indices))
            else:
                for val in inos_in_cluster_contacts:
                    if val > 0:
                        bound_sizes_list.append(1)
                    else:
                        unbound_sizes_list.append(1)

    # compute histograms
    bound_hist = numpy.bincount(numpy.array(bound_sizes_list))
    unbound_hist = numpy.bincount(numpy.array(unbound_sizes_list))
   
    print bound_hist
    print unbound_hist
 
    with open(iso + 'sys' + str(sys) + tag + 'bound_hist.txt', 'w') as results_file:
        for size in range(0, bound_hist.size):
            frequency = bound_hist[size] / float(bound_hist.sum()) 
            results_file.write("%d  %0.2f\n" % (size, frequency))

    with open(iso + 'sys' + str(sys) + tag + 'unbound_hist.txt', 'w') as results_file:
        for size in range(0, unbound_hist.size):
            frequency = unbound_hist[size] / float(unbound_hist.sum())
            results_file.write("%d  %0.2f\n" % (size, frequency))

        
