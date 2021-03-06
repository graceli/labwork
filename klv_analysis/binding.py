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

class BoundState:
    # Starting state
    Start = 1
    # Moved from one state into another. Transition may be irreversible
    Transition = 2

class TransitionError(Exception):
    def __init__(self, value):
        self.value = value
    
    def __str__(self):
        return repr(self.value)


def _binding_constant(polarMatrix, nonpolarMatrix, inositol_concentration):
    total_binding_per_inositol = polarMatrix[:,1:] + nonpolarMatrix[:,1:]
    print "polar matrix nonzero counts", numpy.count_nonzero(numpy.sum(polarMatrix[:, 1:], axis=1))
    print "nonpolar matrix nonzero counts", numpy.count_nonzero(numpy.sum(nonpolarMatrix[:, 1:], axis=1))

    print total_binding_per_inositol 
    total_binding_all_inositols = numpy.sum(total_binding_per_inositol, axis=1)

    print total_binding_all_inositols
 
    bound = numpy.count_nonzero(total_binding_all_inositols)
    unbound = total_binding_all_inositols.size - bound
    print bound, unbound
    binding_constant = float(unbound) * inositol_concentration / bound 
    
    return binding_constant

def _num_binding_events_beta(iso, sys, contacts_ts):
    print contacts_ts[:, 1:]
 
    total_contacts_ts = numpy.sum(contacts_ts[:,1:], axis=1)
    # numpy.savetxt(iso + '_' + str(sys) + '.txt', total_contacts_ts, fmt='%d')

    print total_contacts_ts.shape
    print total_contacts_ts

    num_binding_events = 0
    i = 0
    while i < len(total_contacts_ts):
        while i < len(total_contacts_ts) and total_contacts_ts[i] > 0:
            i += 1
            continue

        if i < len(total_contacts_ts): 
            num_binding_events += 1
        
        while i < len(total_contacts_ts) and total_contacts_ts[i] == 0:
            i += 1
            continue
    
    return num_binding_events


def _num_binding_events(iso, sys, contacts_ts):
    # total_contacts_ts = numpy.sum(contacts_ts[:,1:], axis=1)
    # numpy.savetxt(iso + '_' + str(sys) + '.txt', total_contacts_ts, fmt='%d')

    nrows, ncols = contacts_ts.shape

    print nrows, ncols
    total_num_binding_events = 0
    
    for col in range(0, ncols):
        num_events = 0
        i = 0
        while i < nrows:
            while i < nrows and contacts_ts[i, col] > 0:
                i += 1
                continue

            if i < nrows: 
                num_events += 1

            while i < nrows and contacts_ts[i, col] == 0:
                i += 1
                continue

        total_num_binding_events += num_events

    return total_num_binding_events


def _num_binding_events_state_machine(contacts_ts):
    print contacts_ts.shape
    nrows,ncols = contacts_ts.shape
    total_num_binding_events = 0

    print nrows, ncols

    for col in range(0, ncols):
        print contacts_ts[:, col]
        total_num_binding_events += _count_binding_events_state_machine(contacts_ts[:, col])

    return total_num_binding_events


def _count_binding_events_state_machine(contacts_ts_raw):
    bool_array = contacts_ts_raw >= 1
    total_contacts_ts = bool_array.astype(int)
    
    if len(total_contacts_ts) == 1:
        raise TransitionError("Size of array is one.  No transitions are possible.")

    num_binding_events = 0
    i = 1
    state = BoundState.Start
    while i < len(total_contacts_ts):
        while i < len(total_contacts_ts) and (total_contacts_ts[i] == total_contacts_ts[i-1]):
            i += 1
            continue

        # Hit the end of the time series
        if i == len(total_contacts_ts):
            break

        if state == BoundState.Start:
            state = BoundState.Transition
        else:
            state = BoundState.Start
            num_binding_events += 1

        i += 1

    return num_binding_events


# This function estimates the number of binding events if provided a timeseries
def beta_binding_event_estimate(h5file, inositol_ratio, inositol_concentration, system_indices=[]):
    assert len(system_indices) > 0, "List of system indices should be non-empty."
    
    isomerList = ["scyllo", "chiro"]
    polar_path = "/polar"
    nonpolar_path = "/nonpolar_residue"
    csv_header = ["isomer", "sys_idx", "binding_constant", "inos_conc"]

    writer = csv.writer(open(inositol_ratio + '_binding_events.csv', 'wb'), delimiter=' ')
    writer.writerow(csv_header)

    for iso in isomerList:
        for sys in system_indices[iso]:
            nonpolar_file = os.path.join(nonpolar_path, "%(iso)s_t%(sys)d_per_inositol_contacts.dat" % vars())
            polar_file = os.path.join(polar_path, "%(iso)s_t%(sys)d_inos_total.dat" % vars())
            print polar_file, nonpolar_file

            nonpolar_matrix = myh5.getTableAsMatrix(h5file, nonpolar_file, dtype=numpy.float64)
            polar_matrix = myh5.getTableAsMatrix(h5file, polar_file, dtype=numpy.float64)
            writer.writerow([iso, sys, _num_binding_events_state_machine(nonpolar_matrix[:,1:] + polar_matrix[:,1:]), inositol_concentration])


# This function estimates the number of binding events if provided a timeseries
def disordered_binding_event_estimate(h5file, inositol_ratio, inositol_concentration, system_indices=[]):
    assert len(system_indices) > 0, "List of system indices should be non-empty."
    
    polarName = {'4to2' : 'inos_total.dat', '15to4' : 'whole_nosol_0-200ns_inos_total.dat', '45to4' : 'whole_nosol_0-200ns_inos_total.dat'}
    nonpolarName = {'4to2' : 'per_inositol_contacts.dat', '15to4' : 'whole_nosol_0-200ns_per_inositol_contacts.dat', '45to4' : 'whole_nosol_0-200_per_inositol_contacts.dat'}

    isomerList = ["scyllo", "chiro"]
    polarPath = "/polar"
    nonpolarPath = "/nonpolar_residue"
    csv_header = ["isomer", "sys_idx", "binding_constant", "inos_conc"]
    writer = csv.writer(open(inositol_ratio + '_binding_events.csv', 'wb'), delimiter=' ')
    writer.writerow(csv_header)

    for iso in isomerList:
        data = []
        for sys in system_indices[iso]:
            polarFile = os.path.join(polarPath, "%(iso)s_sys%(sys)s_%(inositol_ratio)s_" % vars() + polarName[inositol_ratio])
            if inositol_ratio == "4to2":
                polarFile = os.path.join(polarPath, "%(iso)s_sys%(sys)s_" % vars() + polarName[inositol_ratio])

            print "analyzing", polarFile

            polarMatrix = myh5.getTableAsMatrix(h5file, polarFile, dtype=numpy.float64)
            print polarMatrix.shape 

            nonpolarFile = os.path.join(nonpolarPath, "%(iso)s_sys%(sys)s_%(inositol_ratio)s_" % vars() + nonpolarName[inositol_ratio])
            if inositol_ratio == "4to2":
                nonpolarFile = os.path.join(nonpolarPath, "%(iso)s_sys%(sys)s_" % vars() + nonpolarName[inositol_ratio])

            print "analyzing", nonpolarFile
            # This is really bad, but for other systems except 4to2 this line was used.  I've changed it 
            if inositol_ratio == "4to2":
                nonpolarMatrix = myh5.getTableAsMatrix(h5file, nonpolarFile, dtype=numpy.float64)[:, :]
            else:
                nonpolarMatrix = myh5.getTableAsMatrix(h5file, nonpolarFile, dtype=numpy.float64)[::2, :]

            print nonpolarMatrix.shape
            
            num_binding_events = _num_binding_events_state_machine(nonpolarMatrix[:,1:] + polarMatrix[:,1:])
            writer.writerow([iso, sys, num_binding_events, inositol_concentration]) 



def compute_disordered_binding_constant(h5file, inositol_ratio, inositol_concentration, system_indices=[]):
    assert len(system_indices) > 0, "List of system indices should be non-empty."
    
    polarName = {'4to2' : 'inos_total.dat', '15to4' : 'whole_nosol_0-200ns_inos_total.dat', '45to4' : 'whole_nosol_0-200_inos_total.dat'}
    nonpolarName = {'4to2' : 'per_inositol_contacts.dat', '15to4' : 'whole_nosol_0-200ns_per_inositol_contacts.dat', '45to4' : 'whole_nosol_0-200_per_inositol_contacts.dat'}

    isomerList = ["scyllo", "chiro"]
    polarPath = "/polar"
    nonpolarPath = "/nonpolar_residue"
    csv_header = ["isomer", "sys_idx", "binding_constant", "inos_conc"]
    writer = csv.writer(open(inositol_ratio + '_binding_constants.csv', 'wb'), delimiter=' ')
    writer.writerow(csv_header)

    for iso in isomerList:
        data = []
        for sys in system_indices[iso]:
            polarFile = os.path.join(polarPath, "%(iso)s_sys%(sys)s_%(inositol_ratio)s_" % vars() + polarName[inositol_ratio])
            if inositol_ratio == "4to2":
                polarFile = os.path.join(polarPath, "klvffae_aggr%(sys)s_%(iso)s_nosol.xtc_" % vars() + polarName[inositol_ratio])

            print "analyzing", polarFile

            polarMatrix = myh5.getTableAsMatrix(h5file, polarFile, dtype=numpy.float64)
            print polarMatrix.shape 

            nonpolarFile = os.path.join(nonpolarPath, "%(iso)s_sys%(sys)s_%(inositol_ratio)s_" % vars() + nonpolarName[inositol_ratio])
            if inositol_ratio == "4to2":
                nonpolarFile = os.path.join(nonpolarPath, "%(iso)s_sys%(sys)s_per_inositol_contacts.dat" % vars())

            print "analyzing", nonpolarFile
            # This is really bad, but for other systems except 4to2 this line was used.  I've changed it 
            if inositol_ratio == "4to2":
                nonpolarMatrix = myh5.getTableAsMatrix(h5file, nonpolarFile, dtype=numpy.float64)[:, :]
            else:
                nonpolarMatrix = myh5.getTableAsMatrix(h5file, nonpolarFile, dtype=numpy.float64)[::2, :]

            print nonpolarMatrix.shape

            binding_constant = _binding_constant(polarMatrix, nonpolarMatrix, inositol_concentration)             
            writer.writerow([iso, sys, binding_constant, inositol_concentration]) 


def compute_beta_binding_constant(h5file, inositol_ratio, inositol_concentration, system_indices=[]):
    assert len(system_indices) > 0, "List of system indices should be non-empty."
    
    isomerList = ["scyllo", "chiro"]
    polar_path = "/polar"
    nonpolar_path = "/nonpolar_residue"
    csv_header = ["isomer", "sys_idx", "binding_constant", "inos_conc"]

    writer = csv.writer(open(inositol_ratio + '_binding_constants.csv', 'wb'), delimiter=' ')
    writer.writerow(csv_header)

    for iso in isomerList:
        for sys in system_indices[iso]:
            nonpolar_file = os.path.join(nonpolar_path, "%(iso)s_t%(sys)d_per_inositol_contacts.dat" % vars())
            polar_file = os.path.join(polar_path, "%(iso)s_t%(sys)d_inos_total.dat" % vars())
            print polar_file, nonpolar_file

            nonpolar_matrix = myh5.getTableAsMatrix(h5file, nonpolar_file, dtype=numpy.float64)
            polar_matrix = myh5.getTableAsMatrix(h5file, polar_file, dtype=numpy.float64)

            print polar_matrix.shape, nonpolar_matrix.shape
            
            binding_constant = _binding_constant(polar_matrix, nonpolar_matrix, inositol_concentration)
            writer.writerow([iso, sys, binding_constant, inositol_concentration])


def compute_beta_low_molar_binding_constant(h5file, inositol_ratio, inositol_concentration):
    # assert len(system_indices) > 0, "List of system indices should be non-empty."
    
    isomerList = ["scyllo", "chiro"]
    polar_path = "/polar"
    nonpolar_path = "/nonpolar_revision"
    csv_header = ["isomer", "sys_idx", "binding_constant", "inos_conc"]

    writer = csv.writer(open('beta_' + inositol_ratio + '_binding_constants.csv', 'wb'), delimiter=' ')
    writer.writerow(csv_header)

    for iso in isomerList:
        for sys in range(0, 3):
            for i in range(1, 6):
                nonpolar_file = os.path.join(nonpolar_path, "%(iso)s_sys%(sys)d_t%(i)d_per_inositol_contacts.dat" % vars())
                polar_file = os.path.join(polar_path, "%(iso)s_sys%(sys)d_t%(i)d_inos_total.dat" % vars())

                print polar_file, nonpolar_file

                nonpolar_matrix = myh5.getTableAsMatrix(h5file, nonpolar_file, dtype=numpy.float64)
                polar_matrix = myh5.getTableAsMatrix(h5file, polar_file, dtype=numpy.float64)

                print polar_matrix.shape, nonpolar_matrix.shape
            
                binding_constant = _binding_constant(polar_matrix, nonpolar_matrix, inositol_concentration)
                writer.writerow([iso, sys, binding_constant, inositol_concentration])


def compute_beta_low_molar_binding_events(h5file, inositol_ratio, inositol_concentration):
    # assert len(system_indices) > 0, "List of system indices should be non-empty."
    
    isomerList = ["scyllo", "chiro"]
    polar_path = "/polar"
    nonpolar_path = "/nonpolar_revision"
    csv_header = ["isomer", "sys_idx", "binding_events", "inos_conc"]

    writer = csv.writer(open('beta_' + inositol_ratio + '_binding_events.csv', 'wb'), delimiter=' ')
    writer.writerow(csv_header)

    for iso in isomerList:
        for sys in range(0, 3):
            for i in range(1, 6):
                nonpolar_file = os.path.join(nonpolar_path, "%(iso)s_sys%(sys)d_t%(i)d_per_inositol_contacts.dat" % vars())
                polar_file = os.path.join(polar_path, "%(iso)s_sys%(sys)d_t%(i)d_inos_total.dat" % vars())

                print polar_file, nonpolar_file

                nonpolar_matrix = myh5.getTableAsMatrix(h5file, nonpolar_file, dtype=numpy.float64)
                polar_matrix = myh5.getTableAsMatrix(h5file, polar_file, dtype=numpy.float64)

                print polar_matrix.shape, nonpolar_matrix.shape
            
                num_binding_events = _num_binding_events_state_machine(polar_matrix[:, 1:] + nonpolar_matrix[:, 1:])
                writer.writerow([iso, sys, num_binding_events, inositol_concentration])

    
def monomer_2to1_binding_events_estimate(h5file, inositol_concentration):
    isomer =  ["scyllo", "chiro"]
   
    csv_header = ["isomer", "sys_idx", "binding_constant", "inos_conc"]
    writer = csv.writer(open('monomer_2to1_binding_events.csv', 'wb'), delimiter=' ')
    writer.writerow(csv_header)

    for iso in isomer:
        for i in range(0, 6):
            polar_matrix = myh5.getTableAsMatrix(h5file, os.path.join('/nonpolar_residue', '%(iso)s_sys%(i)d_mon_2to1_per_inositol_contacts.dat' % vars()), dtype=numpy.float64)
            nonpolar_matrix = myh5.getTableAsMatrix(h5file, os.path.join('/polar', '%(iso)s_sys%(i)d_mon_2to1_inos_total.dat' % vars()), dtype=numpy.float64)

            print polar_matrix.shape
            print nonpolar_matrix.shape
            num_binding_events = _num_binding_events_state_machine(nonpolar_matrix[:, 1:] + polar_matrix[:, 1:])
            writer.writerow([iso, i, num_binding_events, inositol_concentration])

    
def compute_monomer_2to1_binding_constants(h5file, inositol_concentration):
    isomer =  ["scyllo", "chiro"]
   
    csv_header = ["isomer", "sys_idx", "binding_constant", "inos_conc"]
    writer = csv.writer(open('monomer_2to1_binding_constants.csv', 'wb'), delimiter=' ')
    writer.writerow(csv_header)

    for iso in isomer:
        for i in range(0, 6):
            polar_matrix = myh5.getTableAsMatrix(h5file, os.path.join('/nonpolar_residue', '%(iso)s_sys%(i)d_mon_2to1_per_inositol_contacts.dat' % vars()), dtype=numpy.float64)
            nonpolar_matrix = myh5.getTableAsMatrix(h5file, os.path.join('/polar', '%(iso)s_sys%(i)d_mon_2to1_inos_total.dat' % vars()), dtype=numpy.float64)

            print polar_matrix.shape
            print nonpolar_matrix.shape

            binding_constant = _binding_constant(polar_matrix, nonpolar_matrix, inositol_concentration)
            writer.writerow([iso, i, binding_constant, inositol_concentration])

def compute_monomer_15to1_binding_constant(h5file, inositol_ratio, inositol_concentration):
    writer = csv.writer(open('monomer_' + inositol_ratio + '_binding_constants.csv', 'wb'), delimiter=' ')
    csv_header = ["isomer", "inositol_ratio", "binding_constant", "inos_conc"]
    writer.writerow(csv_header)
    for isomer in ["scyllo", "chiro"]:
        for k in range(1, 6):
            polar_big_matrix = None
            nonpolar_big_matrix = None
            from_idx = (k-1)*100 + 1
            to_idx = k*100 + 1
            
            print "Computing run_set", k, "with systems from", from_idx, "to", to_idx-1
            
            for i in range(from_idx, to_idx):
                polar_matrix = myh5.getTableAsMatrix(h5file, '/polar/%(isomer)s_sys%(i)d_inos_total.dat' % vars(), dtype=numpy.float64)
                nonpolar_matrix = myh5.getTableAsMatrix(h5file, '/nonpolar/%(isomer)s_sys%(i)d_per_inositol_contacts.dat' % vars(), dtype=numpy.float64)
               
                if polar_matrix is not None and nonpolar_matrix is not None:
                    if polar_big_matrix is None and nonpolar_big_matrix is None: 
                        polar_big_matrix = polar_matrix
                        nonpolar_big_matrix = nonpolar_matrix
                    else:
                        polar_big_matrix = numpy.concatenate((polar_big_matrix, polar_matrix))
                        nonpolar_big_matrix = numpy.concatenate((nonpolar_big_matrix, nonpolar_matrix))
                else:
                    print "data files (polar and nonpolar) for system", i, "was not found"
            
            binding_constant = _binding_constant(polar_big_matrix, nonpolar_big_matrix, inositol_concentration)
            writer.writerow([isomer, inositol_ratio, binding_constant, inositol_concentration])     
        # binding_constant = _binding_constant(polar_matrix, nonpolar_matrix, inositol_concentration)

def monomer_15to1_binding_events(h5file, inositol_concentration):
    writer = csv.writer(open('monomer_15to1_binding_events.csv', 'wb'), delimiter=' ')
    csv_header = ["isomer", "inositol_ratio", "binding_constant", "inos_conc"]
    writer.writerow(csv_header)
    for isomer in ["scyllo", "chiro"]:
        for k in range(1, 6):
            polar_big_matrix = None
            nonpolar_big_matrix = None
            from_idx = (k-1)*100 + 1
            to_idx = k*100 + 1
            
            print "Computing run_set", k, "with systems from", from_idx, "to", to_idx-1
            
            for i in range(from_idx, to_idx):
                polar_matrix = myh5.getTableAsMatrix(h5file, '/polar/%(isomer)s_sys%(i)d_inos_total.dat' % vars(), dtype=numpy.float64)
                nonpolar_matrix = myh5.getTableAsMatrix(h5file, '/nonpolar/%(isomer)s_sys%(i)d_per_inositol_contacts.dat' % vars(), dtype=numpy.float64)
               
                if polar_matrix is not None and nonpolar_matrix is not None:
                    if polar_big_matrix is None and nonpolar_big_matrix is None: 
                        polar_big_matrix = polar_matrix
                        nonpolar_big_matrix = nonpolar_matrix
                    else:
                        polar_big_matrix = numpy.concatenate((polar_big_matrix, polar_matrix))
                        nonpolar_big_matrix = numpy.concatenate((nonpolar_big_matrix, nonpolar_matrix))
                else:
                    print "data files (polar and nonpolar) for system", i, "was not found"
            
            num_binding_events = _num_binding_events_state_machine(nonpolar_big_matrix[:, 1:] + polar_big_matrix[:, 1:])
            writer.writerow([isomer, "15to1", num_binding_events, inositol_concentration])     


def intersection_mon_low_molar_ratio(h5file, ratio):
    for iso in ['scyllo', 'chiro']:
        for run_set in range(0, 6):
            polar_matrix = myh5.getTableAsMatrix(h5file, os.path.join('/polar', '%(iso)s_sys%(run_set)d_mon_2to1_inos_total.dat' % vars()), dtype=numpy.float64)
            print polar_matrix
            print polar_matrix.shape
            nonpolar_matrix = myh5.getTableAsMatrix(h5file, os.path.join('/nonpolar_residue', '%(iso)s_sys%(run_set)d_mon_2to1_per_inositol_contacts.dat' % vars()), dtype=numpy.float64)
            print nonpolar_matrix
            print nonpolar_matrix.shape

            assert polar_matrix.shape == nonpolar_matrix.shape, "the two matrices are expected to have the same dimensions"

            nrows, ncols = polar_matrix.shape
            counts = [{'polar_only':0, 'nonpolar_only':0, 'polar_nonpolar':0}, {'polar_only':0, 'nonpolar_only':0, 'polar_nonpolar':0}]
            for i in range(0, nrows): 
                for j in range(1, ncols):
                    if polar_matrix[i][j] and nonpolar_matrix[i][j]:
                        counts[0]['polar_nonpolar'] += 1
                    elif polar_matrix[i][j]:
                        counts[0]['polar_only'] += 1
                    elif nonpolar_matrix[i][j]:
                        counts[0]['nonpolar_only'] += 1

            # class csv.DictWriter(csvfile, fieldnames[, restval=''[, extrasaction='raise'[, dialect='excel'[, *args, **kwds]]]])
            total = counts[0]['polar_only'] + counts[0]['nonpolar_only'] + counts[0]['polar_nonpolar']
            fraction = {'polar_only': float(counts[0]['polar_only']) / total, 'nonpolar_only' : float(counts[0]['nonpolar_only']) / total, 'polar_nonpolar' : float(counts[0]['polar_nonpolar']) / total}

            print counts[0]

            writer = csv.DictWriter(open('%(iso)s_sys%(run_set)d_mon_2to1_intersection.csv' % vars(), 'wb'), counts[0].keys())
            writer.writeheader()
            writer.writerow(counts[0])
            writer.writerow(fraction)


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
    polarName = {'4to2' : 'inos_total.dat', '15to4' : 'whole_nosol_0-200ns_inos_total.dat', '45to4' : 'whole_nosol_0-200_inos_total.dat'}
    nonpolarName = {'4to2' : 'per_inositol_contacts.dat', '15to4' : 'whole_nosol_0-200ns_per_inositol_contacts.dat', '45to4' : 'whole_nosol_0-200_per_inositol_contacts.dat'}
   # klvffae_aggr0_chiro_nosol.xtc_inos_total.dat 
    isomerList = ["scyllo", "chiro"]
    polarPath = "/polar"
    nonpolarPath = "/nonpolar_residue"
    dataList = [['isomer', 'system#', 'polar_only', 'polar_and_nonpolar', 'nonpolar_only', 'total']]
    
    resultsWriter = csv.writer(open(ratio + '_intersection.txt', 'wb'), delimiter=' ')
    
    print system_indices
    
    for iso in isomerList:
        data = []
        for sys in system_indices[iso]:
            if ratio == "4to2":
                polarFile = os.path.join(polarPath, "klvffae_aggr%(sys)s_%(iso)s_nosol.xtc_" % vars() + polarName[ratio])
            else:
                polarFile = os.path.join(polarPath, "%(iso)s_sys%(sys)s_%(ratio)s_" % vars() + polarName[ratio])

            print "analyzing", polarFile
            polarMatrix = myh5.getTableAsMatrix(h5file, polarFile, dtype=numpy.float64)
            print polarMatrix
            
            if ratio == "4to2":
                nonpolarFile = os.path.join(nonpolarPath, "%(iso)s_sys%(sys)s_" % vars() + nonpolarName[ratio])
            else:
                nonpolarFile = os.path.join(nonpolarPath, "%(iso)s_sys%(sys)s_%(ratio)s_" % vars() + nonpolarName[ratio])
            print "analyzing", nonpolarFile

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
        std = numpy.std(numpy.array(data), axis=0) / len(system_indices)
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

def intersection_beta_low_molar(h5file):
    for iso in ["scyllo", "chiro"]:
        for sys in range(0, 3):
            for i in range(1, 6):
                nonpolar_file = os.path.join('/nonpolar_revision', "%(iso)s_sys%(sys)d_t%(i)d_per_inositol_contacts.dat" % vars())
                polar_file = os.path.join('/polar', "%(iso)s_sys%(sys)d_t%(i)d_inos_total.dat" % vars())
                _intersection(h5file, polar_file, nonpolar_file, "%(iso)s_beta_low_molar_ratio_sys%(sys)d_t%(i)d" % vars())
                
def _intersection(h5file, polar_file, nonpolar_file, tag):
    nonpolar_matrix = myh5.getTableAsMatrix(h5file, nonpolar_file, dtype=numpy.float64)
    polar_matrix = myh5.getTableAsMatrix(h5file, polar_file, dtype=numpy.float64)

    counts = {'polar_only' : 0, 'nonpolar_only' : 0, 'polar_nonpolar' : 0}

    if polar_matrix is not None and nonpolar_matrix is not None:
        print polar_matrix.shape
        print nonpolar_matrix.shape

        assert polar_matrix.shape == nonpolar_matrix.shape, "the two matrices are expected to have the same dimensions"

        nrows, ncols = polar_matrix.shape
        for i in range(1, nrows):
            for j in range(1, ncols):
                if polar_matrix[i][j] and nonpolar_matrix[i][j]:
                    counts['polar_nonpolar'] += 1
                elif polar_matrix[i][j]:
                    counts['polar_only'] += 1
                elif nonpolar_matrix[i][j]:
                    counts['nonpolar_only'] += 1

    # normalize
    total = counts['polar_nonpolar'] + counts['polar_only'] + counts['nonpolar_only']
    if total != 0:
        counts['polar_nonpolar'] = counts['polar_nonpolar'] / float(total)
        counts['polar_only'] = counts['polar_only'] / float(total)
        counts['nonpolar_only'] = counts['nonpolar_only'] / float(total)


    writer = csv.DictWriter(open('%(tag)s_intersection.csv' % vars(), 'wb'), counts.keys())
    writer.writeheader()
    write_header = True

    # class csv.DictWriter(csvfile, fieldnames[, restval=''[, extrasaction='raise'[, dialect='excel'[, *args, **kwds]]]])
    print sys, counts
    writer.writerow(counts)
               
def intersection_beta(h5file, tag):
    isomerList = ["scyllo", "chiro"]
    polar_path = "/polar"
    nonpolar_path = "/nonpolar_residue"
   
    write_header = False
    for iso in isomerList:
        for sys in range(0, 6):
            nonpolar_file = os.path.join(nonpolar_path, "%(iso)s_t%(sys)d_per_inositol_contacts.dat" % vars())
            polar_file = os.path.join(polar_path, "%(iso)s_t%(sys)d_inos_total.dat" % vars())

            nonpolar_matrix = myh5.getTableAsMatrix(h5file, nonpolar_file, dtype=numpy.float64)
            polar_matrix = myh5.getTableAsMatrix(h5file, polar_file, dtype=numpy.float64)

            counts = {'polar_only' : 0, 'nonpolar_only' : 0, 'polar_nonpolar' : 0}

            if polar_matrix is not None and nonpolar_matrix is not None:
                print polar_matrix.shape
                print nonpolar_matrix.shape

                assert polar_matrix.shape == nonpolar_matrix.shape, "the two matrices are expected to have the same dimensions"

                nrows, ncols = polar_matrix.shape
                for i in range(1, nrows):
                    for j in range(1, ncols):
                        if polar_matrix[i][j] and nonpolar_matrix[i][j]:
                            counts['polar_nonpolar'] += 1
                        elif polar_matrix[i][j]:
                            counts['polar_only'] += 1
                        elif nonpolar_matrix[i][j]:
                            counts['nonpolar_only'] += 1

            # normalize
            total = counts['polar_nonpolar'] + counts['polar_only'] + counts['nonpolar_only']
            if total != 0:
                counts['polar_nonpolar'] = counts['polar_nonpolar'] / float(total)
                counts['polar_only'] = counts['polar_only'] / float(total)
                counts['nonpolar_only'] = counts['nonpolar_only'] / float(total)


            writer = csv.DictWriter(open('%(iso)s_%(tag)s_intersection%(sys)d.csv' % vars(), 'wb'), counts.keys())
            if write_header is False:
                writer.writeheader()
            else:
                write_header = True

            # class csv.DictWriter(csvfile, fieldnames[, restval=''[, extrasaction='raise'[, dialect='excel'[, *args, **kwds]]]])
            print sys, counts
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


def nonpolar_residue_monomer_2to1(tag):
    for isomer in ['scyllo', 'chiro']:
        data_list = []
        for i in range(0, 6):
            # chiro_run_set4_per_residue_contacts.dat
            file = "%(isomer)s_run_set%(i)d_per_residue_contacts.dat" % vars()
            print "analyzing ", file

            # read in the file
            data = numpy.genfromtxt(file, comments="#", dtype='float')
            nrows,ncols = data.shape
            print nrows

            # sum over rows
            time_avg = numpy.average(data[:,1:], axis=0)
            print time_avg
            print time_avg.shape
            # time_avg.shape = (time_avg.size, ) 
            # print time_avg.shape
            # sum_over_peptides = numpy.sum(time_avg, axis=1)

            data_list.append(time_avg)

        # save results to flat files
        nparray = numpy.array(data_list)

        # dump the list of results for each system
        numpy.savetxt('%(isomer)s_%(tag)s_nonpolar_residue_contact.txt' % vars(), nparray, fmt='%0.8f')

        # average over all the systems; each system is a row in nparray
        average = numpy.average(nparray, axis=0)
        std = numpy.std(nparray, axis=0) / math.sqrt(5)

        #save the normalized average and std
        numpy.savetxt('%(isomer)s_%(tag)s_nonpolar_residue_contact_avg_std.txt' % vars(), [average, std], fmt='%0.8f')


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

            print table_path

            data = myh5.getTableAsMatrix(h5file, table_path, dtype=numpy.float64)
            
            # print data
            
            sum_over_time = numpy.average(data[20000:N_datapoints, 1:], axis = 0)
            
            # print sum_over_time
            
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


def nonpolar_residue_beta_low_molar(h5file, system_indices=[]):
    assert len(system_indices) > 0, "The list of system_indices should not be empty"

    for isomer in ["scyllo", "chiro"]:
        data_list = []
        for s in system_indices:
            nonpolar_residue_large = None
            for i in range(1,6):
                nonpolar_residue_path = "/nonpolar_revision/%(isomer)s_sys%(s)d_t%(i)d_per_residue_contacts.dat" % vars()
                print "analyzing ", nonpolar_residue_path

                # read in the file
                data = myh5.getTableAsMatrix(h5file, nonpolar_residue_path, dtype=numpy.float64)
                print data
 
                if nonpolar_residue_large is None:
                    print "in here"
                    nonpolar_residue_large = data
                else:
                    nonpolar_residue_large = numpy.concatenate((nonpolar_residue_large, data))
                
                # data = numpy.genfromtxt(file, comments="#", dtype='float')
                nrows,ncols = nonpolar_residue_large.shape
                print nrows, ncols

            # sum over rows
            time_avg = numpy.average(nonpolar_residue_large[:,1:], axis=0)
            print time_avg.shape
            time_avg.shape = (time_avg.size / 16, 16) 
            print time_avg.shape
            sum_over_peptides = numpy.sum(time_avg, axis = 1)

            data_list.append(sum_over_peptides)

        # save results to flat files
        nparray = numpy.array(data_list)

        # dump the list of results for each system
        numpy.savetxt('%(isomer)s_low_molar_nonpolar_residue_contact.txt' % vars(), nparray, fmt='%0.8f')

        # average over all the systems; each system is a row in nparray
        average = numpy.average(nparray, axis=0) / 16 
        std = numpy.std(nparray, axis=0) / 16 / math.sqrt(len(system_indices))

        #save the normalized average and std
        numpy.savetxt('%(isomer)s_low_molar_nonpolar_residue_contact_avg_std.txt' % vars(), [average, std], fmt='%0.8f')


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
