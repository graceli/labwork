#this file defines all the different types of table rows 
#each corresponds to one type of analysis
from tables import *

SASTable = {
	'temp' : Int32Col(dflt=0, pos=0),
	'replicanum' : Int32Col(dflt=0, pos=1),
	'seqnum' : Int32Col(dflt=0, pos=2),
	'time' : Int32Col(dflt=0, pos=3),
	'hydrophobic' : Float32Col(dflt=0.0, pos=4),
	'hydrophilic' : Float32Col(dflt=0.0, pos=5),
	'total' : Float32Col(dflt=0.0, pos=6)
}

RGTable = {
	'temp' : Int32Col(dflt=0, pos=0),
	'replicanum' : Int32Col(dflt=0, pos=1),
	'seqnum' : Int32Col(dflt=0, pos=2),
	'time' : Int32Col(dflt=0, pos=3),
	'Rg' : Float32Col(dflt=0.0, pos=4),
	'Rgx' : Float32Col(dflt=0.0, pos=5),
	'Rgy' : Float32Col(dflt=0.0, pos=6),
	'Rgz' : Float32Col(dflt=0.0, pos=7)
}

EETable = {
	'temp' : Int32Col(dflt=0, pos=0),
	'replicanum' : Int32Col(dflt=0, pos=1),
	'seqnum' : Int32Col(dflt=0, pos=2),
	'time' : Float32Col(dflt=0, pos=3),
	'eed' : Float32Col(dflt=0.0, pos=4)
}

RamaTable = {
	'temp' : Int32Col(dflt=0, pos=0),
	'replicanum' : Int32Col(dflt=0, pos=1),
	'seqnum' : Int32Col(dflt=0, pos=2),
	#'time' : Float32Col(dflt=0, pos=3),
	'phi' : Float32Col(dflt=0.0, pos=4),
	'psi' : Float32Col(dflt=0.0, pos=5),
	'residue' : StringCol(16,pos=6)
}

QTable = {
	'temp' : Int32Col(dflt=0, pos=0),
	'replicanum' :  Int32Col(dflt=0, pos=1),
	'seqnum' : Int32Col(dflt=0, pos=2),
	'framenative' : Int32Col(dflt=0, pos=3),
	'total' : Int32Col(dflt=0, pos=4),
	'nativetotal' : Int32Col(dflt=0, pos=5),
	'Q' : Float32Col(dflt=0, pos=6)
}

ContactMapTable = {
	'temp' : Int32Col(dflt=0, pos=0),
	'replicanum' : Int32Col(dflt=0,pos=1),
	'seqnum' : Int32Col(dflt=0, pos=2),
	'res1' : Int32Col(dflt=0, pos=3), 'res2' : Int32Col(dflt=0, pos=4), 'res3' : Int32Col(dflt=0, pos=5),
	'res4' : Int32Col(dflt=0, pos=6), 'res5' : Int32Col(dflt=0, pos=7), 'res6' : Int32Col(dflt=0, pos=8),
	'res7' : Int32Col(dflt=0, pos=9), 'res8' : Int32Col(dflt=0, pos=10), 'res9' : Int32Col(dflt=0, pos=11),
	'res10' : Int32Col(dflt=0, pos=12), 'res11' : Int32Col(dflt=0, pos=13), 'res12' : Int32Col(dflt=0, pos=14),
	'res13' : Int32Col(dflt=0, pos=15), 'res14' : Int32Col(dflt=0, pos=16), 'res15' : Int32Col(dflt=0, pos=17),
	'res16' : Int32Col(dflt=0, pos=18), 'res17' : Int32Col(dflt=0, pos=19), 'res18' : Int32Col(dflt=0, pos=20),
	'res19' : Int32Col(dflt=0, pos=21), 'res20' : Int32Col(dflt=0, pos=22), 'res21' : Int32Col(dflt=0, pos=23), 
	'res22' : Int32Col(dflt=0, pos=24), 'res23' : Int32Col(dflt=0, pos=25), 'res24' : Int32Col(dflt=0, pos=26), 
	'res25' : Int32Col(dflt=0, pos=27), 'res26' : Int32Col(dflt=0, pos=28), 'res27' : Int32Col(dflt=0, pos=29), 
	'res28' : Int32Col(dflt=0, pos=30), 'res29' : Int32Col(dflt=0, pos=31), 'res30' : Int32Col(dflt=0, pos=32), 
	'res31' : Int32Col(dflt=0, pos=33), 'res32' : Int32Col(dflt=0, pos=34), 'res33' : Int32Col(dflt=0, pos=35), 
	'res34' : Int32Col(dflt=0, pos=36), 'res35' : Int32Col(dflt=0, pos=37), 'res36' : Int32Col(dflt=0, pos=38), 
	'res37' : Int32Col(dflt=0, pos=39), 'res38' : Int32Col(dflt=0, pos=40), 'res39' : Int32Col(dflt=0, pos=41), 
	'res40' : Int32Col(dflt=0, pos=42), 'res41' : Int32Col(dflt=0, pos=43), 'res42' : Int32Col(dflt=0, pos=44), 
	'res43' : Int32Col(dflt=0, pos=45), 'res44' : Int32Col(dflt=0, pos=46), 'res45' : Int32Col(dflt=0, pos=47), 
	'res46' : Int32Col(dflt=0, pos=48), 'res47' : Int32Col(dflt=0, pos=49), 'res48' : Int32Col(dflt=0, pos=50), 
	'res49' : Int32Col(dflt=0, pos=51), 'res50' : Int32Col(dflt=0, pos=51), 'res51' : Int32Col(dflt=0, pos=52), 
	'res52' : Int32Col(dflt=0, pos=53), 'res53' : Int32Col(dflt=0, pos=54), 'res54' : Int32Col(dflt=0, pos=55), 
	'res55' : Int32Col(dflt=0, pos=56), 'res56' : Int32Col(dflt=0, pos=57), 'res57' : Int32Col(dflt=0, pos=58), 
	'res58' : Int32Col(dflt=0, pos=60), 'res59' : Int32Col(dflt=0, pos=61)
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

