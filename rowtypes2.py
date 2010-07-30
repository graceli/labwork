#this file defines all the different types of table rows 
#each corresponds to one type of analysis
from tables import *

#SASTable = {
#	'time' : Int32Col(),
#	'replicanum' : Int32Col(),
#	'seqnum' : Int32Col(),
#	'temp' : Int32Col(),
#	'hydrophobic' : Float32Col(),
#	'hydrophilic' : Float32Col(),
#	'total' : Float32Col()
#}
#
#RGTable = {
#	'time' : Int32Col(),
#	'replicanum' : Int32Col(),
#	'seqnum' : Int32Col(),
#	'temp' : Int32Col(),
#	'Rg' : Float32Col()
#}
#
DefaultTable = {
	'time' : Int32Col(pos=1),
	'replicanum' : Int32Col(pos=2),
	'seqnum' : Int32Col(pos=3),
	'temp' : Int32Col(pos=4),
}

#class DefaultTable(IsDescription):
#	replicanum = Int32Col(pos=1)
#	seqnum = Int32Col(pos=3)
#	time = Int32Col(pos=2)
#	temp = Int32Col(pos=4)
	

if __name__ == "__main__":	
	
	d = Description(DefaultTable)
	print d	
	print d._v_names
