#######################################  basic #######################################  
# iterate over a list
squares = [1, 4, 9 , 16]
sum = 0 
for num in squares:
	sum+=num

print sum


# check for an element in a list
list = ['larry', 'curly', 'moe']
if 'curly' in list:
	print 'yay'
else:
	print 'no'
	

#generate a sequence of numbers
for i in range(100):
	print i


#list slicing
list = ['a', 'b', 'c', 'd']
print list[1:-1]  #['b', 'c']
list[0:2] = 'z'
print list

## Can build up a dict by starting with the the empty dict {}
 ## and storing key/value pairs into the dict like this:
 ## dict[key] = value-for-that-key
 dict = {}
 dict['a'] = 'alpha'
 dict['g'] = 'gamma'
 dict['o'] = 'omega'

 print dict  ## {'a': 'alpha', 'o': 'omega', 'g': 'gamma'}

 print dict['a']     ## Simple lookup, returns 'alpha'
 dict['a'] = 6       ## Put new key/value into dict
 'a' in dict         ## True
 ## print dict['z']                  ## Throws KeyError
 if 'z' in dict: print dict['z']     ## Avoid KeyError
 print dict.get('z')  ## None (instead of KeyError)


for key in dict:
		print key
		
print dict.keys()	## ['a','o','g']
print dict.values() ## ['alpha', 'omega', 'gamma']
print dict.items()  ## [('a', 'alpha'),('o', 'omega'), ('g', 'gamma'')]

del dict['a']
print dict			## {'o': 'omega', 'g': 'gamma'}


#######################################  numpy ####################################### 

from numpy import *

a = array([1,2,3])
b = array((10,11,12))

>>> print a+b
array([11,13,15])

>>> print a.dtype
dtype('<i4')

a = array([1,2,3], dtype=float)

t = arange(0, 2*pi, 0.1)
sinvalues = sin(t)

# slicing
>>>t[:]              # get all t-values
array([ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1. ,
        1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2. ,  2.1,
        2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,  3. ,  3.1,  3.2,
        3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,  4. ,  4.1,  4.2,  4.3,
        4.4,  4.5,  4.6,  4.7,  4.8,  4.9,  5. ,  5.1,  5.2,  5.3,  5.4,
        5.5,  5.6,  5.7,  5.8,  5.9,  6. ,  6.1,  6.2])

>>>t[2:4]            # get sub-array with the elements at the indexes 2,3
array([  0.2,  0.3,  0.4 ])


>>>t[0:6:2]          # every even-indexed value up to but excluding 6
array([ 0. ,  0.2,  0.4])


# 2D arrays
>>> b = arange(12).reshape(3,4)
>>> b
array([[ 0,  1,  2,  3],
       [ 4,  5,  6,  7],
       [ 8,  9, 10, 11]])
>>>
>>> b.sum(axis=0)                            # sum of each column
array([12, 15, 18, 21])
>>>
>>> b.min(axis=1)                            # min of each row
array([0, 4, 8])
>>>
>>> b.cumsum(axis=1)                         # cumulative sum along the rows
array([[ 0,  1,  3,  6],
       [ 4,  9, 15, 22],
       [ 8, 17, 27, 38]])

#2D slicing
>>> t_mat = t.reshape(t.size/9, 9)
>>> t_mat
array([[ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8],
       [ 0.9,  1. ,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7],
       [ 1.8,  1.9,  2. ,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6],
       [ 2.7,  2.8,  2.9,  3. ,  3.1,  3.2,  3.3,  3.4,  3.5],
       [ 3.6,  3.7,  3.8,  3.9,  4. ,  4.1,  4.2,  4.3,  4.4],
       [ 4.5,  4.6,  4.7,  4.8,  4.9,  5. ,  5.1,  5.2,  5.3],
       [ 5.4,  5.5,  5.6,  5.7,  5.8,  5.9,  6. ,  6.1,  6.2]])
>>> t_mat[0]
array([ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8])

>>> t_mat[:,0:2]
array([[ 0. ,  0.1],
       [ 0.9,  1. ],
       [ 1.8,  1.9],
       [ 2.7,  2.8],
       [ 3.6,  3.7],
       [ 4.5,  4.6],
       [ 5.4,  5.5]])

>>> t_mat[0:3, 0:5] 
array([[ 0. ,  0.1,  0.2,  0.3,  0.4],
       [ 0.9,  1. ,  1.1,  1.2,  1.3],
       [ 1.8,  1.9,  2. ,  2.1,  2.2]])

#######################################  plotting #######################################  
"""
Make a histogram of normally distributed random numbers and plot the
analytic PDF over it
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

mu, sigma = 100, 15
x = mu + sigma * np.random.randn(10000)

fig = plt.figure()
ax = fig.add_subplot(111)

# the histogram of the data
n, bins, patches = ax.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

# hist uses np.histogram under the hood to create 'n' and 'bins'.
# np.histogram returns the bin edges, so there will be 50 probability
# density values in n, 51 bin edges in bins and 50 patches.  To get
# everything lined up, we'll compute the bin centers
bincenters = 0.5*(bins[1:]+bins[:-1])
# add a 'best fit' line for the normal PDF
y = mlab.normpdf( bincenters, mu, sigma)
l = ax.plot(bincenters, y, 'r--', linewidth=1)

ax.set_xlabel('Smarts')
ax.set_ylabel('Probability')
#ax.set_title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
ax.set_xlim(40, 160)
ax.set_ylim(0, 0.03)
ax.grid(True)

plt.show()

####################################### plotting code ########################################


data = randn(100000)
hist(data, 50, normed=1, facecolor='green', alpha=0.75)
xlabel('x')
ylabel('Probability')
grid(True)

########################################## Querying ######################################
file = tables.openFile(analysisfile)
table = file.getNode('/'+group+'/'+tablename)

for T in Tlist:
	readout = table.readWhere('temp==%(T)d' % vars())
	
############################################################################################




#######################################  Putting it all together #######################################

import numpy
import tables
import pylab

table_row_descr = {
	'time':Float32Col(dflt=0.0, pos=0),
 	'rmsd':Float32Col(dflt=0.0, pos=1)
}

group_name = 'protein'
table_name = 'rmsd'
table_path = '/protein/rmsd'
h5_filename = 'analysis.h5'

h5file = tables.openFile(h5_filename, mode="a")				# open or create a hdf5 file
filters = tables.Filters(complevel=8, complib='zlib')		# specify compression 

if not h5file.__contains__('/' + group_name):
	h5file.createGroup(h5file.root, group_name, filters=filters)


if not h5file.__contains__(table_path):
	print "table does not exist; creating table", table_name
	table = h5file.createTable('/' + group_name, table_name, table_row_descr)
else:
	table = h5file.getNode(table_path)

data = numpy.genfromtxt('rmsd.dat', dtype=float)

print data
print data.shape

print "max rmsd", numpy.max(data)
print "min rmsd", numpy.min(data)
print "average rmsd", numpy.average(data[:,1], axis=0)

#append numpy array to table
table.append(data)

readout = table.read()

pylab.plot(readout)


















