import os
import numpy
import plot_and_save2hdf5 as myh5


# For each analysis need some preprocessing before saving into a table
# Each preprocessing step will be different depending on the analysis
# Write on function that will do a different preprocessing depending on the 'analysis' type


def generate_file_name(ratio, isomer, sys_idx, analysis):
    if analysis == "rmsd":
        return "{0}/{1}/{3}/{2}_{3}.xvg".format(ratio, isomer, sys_idx, analysis)

    
def rmsd(**kwargs):
    print kwargs, "in rmsd"
    data = kwargs['data']
    if kwargs['keep_time']:
        return data

    return data[:,1]

def preprocess(function, **kwargs):
    """ Reads in and preprocesses the data file as a numpy array 
        according to its analysis type """

    print kwargs, "preprocess"
    return function(**kwargs)
        
def load_into_tables():
    """docstring for load_into_tables"""
    # Initialize a pytable for saving data into

    function_list = {'rmsd' : rmsd } 

    print "initializing pytables data file"
    
    group_name = 'analysis'

    h5file = myh5.initialize('analysis.h5', group_name)

    root = '/' + group_name

    print "reading in flat files"    
    for ratio in [15, 64]:
        for isomer in ["chiro", "scyllo", "glycerol"]:
            for analysis in ["rmsd"]:
                for sys_idx in range(0, 10):
                    flat_file_path = generate_file_name(ratio, isomer, sys_idx, analysis)
                    print "loading in file at", flat_file_path

                    flat_file_name = flat_file_path.replace('/','_')
                    print "loading in", flat_file_name
                    if os.path.exists(flat_file_path):
                        data_file = numpy.genfromtxt(flat_file_path)
                    else:
                        print flat_file_path, "was not found!"
                    
                    data_cleaned = preprocess(function_list[analysis], data=data_file, keep_time=True)
                    #kwargs={'data': data_file, 'keep_time': True})
                    myh5.save(h5file, data_cleaned, os.path.join(root, os.path.splitext(flat_file_name)[0]))
                
                
def main():
    """docstring for main"""
    load_into_tables()
                
    
if __name__ == '__main__':
    main()
