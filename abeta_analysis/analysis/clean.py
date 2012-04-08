import subprocess
import os
import numpy
import tables

# paste together the columns and remove the extra time columns for 2 to RATIO
# t inos1 inos2 ... inosRATIO

# analysisName.h5:/isomer_ratio_analysis_index
# hbonds.h5
# nonpolar.h5

isomer_list = ["chiro"]
ratio_list = [15]

def process_hbonds_inositol(h5file):
    for isomer in isomer_list:
        for ratio in ratio_list:
            for sys in range(10):
                data_path = os.path.join(str(ratio) + '/' + isomer + '/' + 'hbonds_inositol', str(sys))
                files = [ os.path.join(data_path, str(i) + '.xvg') for i in range(1, ratio+1) ]
                print files
                first_file = True
                column_stack = []
                for f in files:
                    print f
                    data = numpy.genfromtxt(f)
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
                file_name = isomer + '_' + str(ratio) + '_' + 'hbonds_inositol' + '_' + str(sys)
                data_all = numpy.transpose(numpy.vstack(column_stack))
                print data_all
                print "saving data into", h5file.root, file_name
                atom = tables.Atom.from_dtype(data_all.dtype)
                filters = tables.Filters(complib='zlib', complevel = 5)
                h5_carray = h5file.createCArray(h5file.root, file_name, atom, data_all.shape, filters=filters)
                h5_carray[:] = data_all

def process_nonpolar_inositol(h5file):
    """docstring for process_nonpolar_inositol"""
    # read in each of the chains and sum up the matrices
    for isomer in isomer_list:
        for ratio in ratio_list:
            for sys in range(10):
                data_all = None
                for ch in range(5):
                    data_path = os.path.join(str(ratio) + '/' + isomer + '/' + 'nonpolar', '%(sys)s_chain%(ch)s_inositol_np_contact' % vars() + '.dat')
                    a_chain = numpy.genfromtxt(data_path)
                    if ch == 0:
                        data_all = a_chain
                    else:
                        data_all[:, 1:] += a_chain[:,1:]

                # data_all = numpy.transpose(numpy.vstack(column_stack))
                # print data_all
                file_name = isomer + '_' + str(ratio) + '_' + 'nonpolar' + '_' + str(sys)
                print "saving data into", h5file.root, file_name
                atom = tables.Atom.from_dtype(data_all.dtype)
                filters = tables.Filters(complib='zlib', complevel = 5)
                h5_carray = h5file.createCArray(h5file.root, file_name, atom, data_all.shape, filters=filters)
                h5_carray[:] = data_all
        
def main():
    h5file = tables.openFile("analysis.h5", 'w')
    print "reading"
    process_hbonds_inositol(h5file)
    process_nonpolar_inositol(h5file)
    h5file.close()

if __name__ == '__main__':
    main()
