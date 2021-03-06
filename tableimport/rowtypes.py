#this file defines all the different types of table rows 
#each corresponds to one type of analysis
from tables import *

SASTable = {
	'temp' : Int32Col(dflt=0, pos=0),
	'replicanum' : Int32Col(dflt=0, pos=1),
	'seqnum' : Int32Col(dflt=0, pos=2),
	'partnum' : Int32Col(dflt=0, pos=3),
	'time' : Int32Col(dflt=0, pos=4),
	'hydrophobic' : Float32Col(dflt=0.0, pos=5),
	'hydrophilic' : Float32Col(dflt=0.0, pos=6),
	'total' : Float32Col(dflt=0.0, pos=7)
}

RGTable = {
	'temp' : Int32Col(dflt=0, pos=0),
	'replicanum' : Int32Col(dflt=0, pos=1),
	'seqnum' : Int32Col(dflt=0, pos=2),
	'partnum' : Int32Col(dflt=0, pos=3),
	'time' : Int32Col(dflt=0, pos=4),
	'Rg' : Float32Col(dflt=0.0, pos=5),
	'Rgx' : Float32Col(dflt=0.0, pos=6),
	'Rgy' : Float32Col(dflt=0.0, pos=7),
	'Rgz' : Float32Col(dflt=0.0, pos=8)
}


DefaultTable = {
	'time' : Int32Col(pos=0),
	'replicanum' : Int32Col(pos=1),
	'seqnum' : Int32Col(pos=2),
	'temp' : Int32Col(pos=3),
}


if __name__ == "__main__":	
	d = Description(DefaultTable)
	print d	
	print d._v_names
	print d._v_types

	r = Description(RGTable)
	print r
	print r._v_names
	print r._v_types

	s = Description(SASTable)
	print s
	print s._v_names
	print s._v_types
