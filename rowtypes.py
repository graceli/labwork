#this file defines all the different types of table rows 
#each corresponds to one type of analysis
from tables import *

SASTable = {
	'time' : Int32Col(),
	'replicanum' : Int32Col(),
	'seqnum' : Int32Col(),
	'temp' : Int32Col(),
	'hydrophobic' : Float32Col(),
	'hydrophilic' : Float32Col(),
	'total' : Float32Col()
}

RGTable = {
	'time' : Int32Col(),
	'replicanum' : Int32Col(),
	'seqnum' : Int32Col(),
	'temp' : Int32Col(),
	'Rg' : Float32Col()
}

DefaultTable = {
	'time' : Int32Col(),
	'replicanum' : Int32Col(),
	'seqnum' : Int32Col(),
	'temp' : Int32Col(),
}