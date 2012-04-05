import numpy
import plot_and_save2hdf5 as myh5


# For each analysis need some preprocessing before saving into a table
# Each preprocessing step will be different depending on the analysis
# Write on function that will do a different preprocessing depending on the 'analysis' type


def generate_file_name(ratio, isomer, sys_idx, analysis):
    if analysis == "rmsd":
        return "{0}/{1}/{2}_{3}.xvg".format(ratio, isomer, sys_idx, analysis)

def rmsd(data, **kwargs):
    if kwargs['keep_time']:
        return data

    return data[:,1]

def preprocess(function, data, **kwargs):
    """ Reads in and preprocesses the data file as a numpy array 
        according to its analysis type """
        
    return function(data, kwargs)
        
def load_into_tables():
    """docstring for load_into_tables"""
    # Initialize a pytable for saving data into
    h5file = myh5.initialize('analysis.h5','/analysis')
    
    for ratio in [15, 64]:
        for isomer in ["chiro", "scyllo", "glycerol"]:
            for analysis in ["rmsd", "rmsf", "chain_hbonds"]
                for sys_idx in range(0, 10):
                    flat_file_name = generate_file_name(ratio, isomer, sys_idx, analysis)
                    data_file = numpy.genfromtxt(flat_filename)
                    data_cleaned = preprocess(analysis, datafile)
                    myh5.save(h5file, data_cleaned, '/' + flat_file_name)
                
                
def main():
    """docstring for main"""
    load_into_tables()
                
    
if __name__ == '__main__':
    main()