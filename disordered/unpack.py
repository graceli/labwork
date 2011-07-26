# Flow for grabbing only replicas for T=300K
# Use a large memory node
# unpack tar
# pick out replica xtcs and edrs at temperatures only at 300K 
#	extract corresponding edr file for xtc
# copy to a new directory in memory
# tar files and copy back down to scratch 
# end result is a bunch of tar files each containing small trajectories at temperature 300K

# read the code for analyse_force_database

# comments from the DR code base
# w_nominal is fixed and stores the start position of each replica
# w refers to the coordinate. w can be a spatial coordinate (which can be the fourth dimension), 
# or w can be lambda, or can be beta if a temperature coordinate is used.

# w_nominal -- the index of a replica. start from 0 to N-1, for N replicas.
# Note that for the case of temperature, if T_i > T_j for i<j, then w_nominal(T_i) = w_nominal(T_j)-1.  eg. 0: 640, 1: 635
# w -- the actual value that is being varied. In the case of temperature, its beta = 1/(R*T). 


# Easily read metadata for a h5 file using the pytables package
# Pytables is not as easy to use ... 
# is there a wrapper tool where I can just dump the column, table information other than ptdump?

