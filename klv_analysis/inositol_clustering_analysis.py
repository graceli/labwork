import csv

# Read in the corresponding file for inositol bound / unbound data

# Parse the csv inositol clusters data
with open('clust_info.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        indices = row[2].split(' ')
        print indices
        # convert the indices to integers using map and a lambda function
        
        # Use the inositol residue indices of the inositols found in the cluster to get the elements from the inositol bound and unbound data
        # Note that I will need to subtract the offset, the starting residue index of inositol molecules in the gro file - eg. 145 for beta, from the indices in the cluster.  The indices in the cluster are inositol residue numbers in the gro file.
        # the indices in the ub/b data file is indexed by the column. ie. 0th column is the 0th inositol corresponding to inositol residue id of 145.
        
        # use the list in the numpy array slicing eg
        # In [3]: a = numpy.array([1,2,3,4])
        # 
        # In [4]: a
        # Out[4]: array([1, 2, 3, 4])
        # 
        # In [5]: a[[1,2]]
        # Out[5]: array([2, 3])
        
        # Sum the resulting array.  If the sum is greater than 0, then this means that at least one inositol in the cluster is bound, and therefore the cluster is bound.
        # Compute the size of the array and bin it for the bound histogram
        # Otherwise, do the same but bin the result in the unbound cluster histogram
    
    # Output the results
    # plot as well?