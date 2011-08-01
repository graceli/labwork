#this file defines all the different types of table rows 
#each corresponds to one type of analysis
from tables import *

# SASTable = {
# 	'temp' : Int64Col(dflt=0, pos=0),
# 	'replicanum' : Int64Col(dflt=0, pos=1),
# 	'seqnum' : Int64Col(dflt=0, pos=2),
# 	'time' : Int64Col(dflt=0, pos=3),
# 	'hydrophobic' : Float64Col(dflt=0.0, pos=4),
# 	'hydrophilic' : Float64Col(dflt=0.0, pos=5),
# 	'total' : Float64Col(dflt=0.0, pos=6)
# }
# 
# RGTable = {
# 	'temp' : Int64Col(dflt=0, pos=0),
# 	'replicanum' : Int64Col(dflt=0, pos=1),
# 	'seqnum' : Int64Col(dflt=0, pos=2),
# 	'time' : Int64Col(dflt=0, pos=3),
# 	'Rg' : Float64Col(dflt=0.0, pos=4),
# 	'Rgx' : Float64Col(dflt=0.0, pos=5),
# 	'Rgy' : Float64Col(dflt=0.0, pos=6),
# 	'Rgz' : Float64Col(dflt=0.0, pos=7)
# }
# 
# EETable = {
# 	'temp' : Int64Col(dflt=0, pos=0),
# 	'replicanum' : Int64Col(dflt=0, pos=1),
# 	'seqnum' : Int64Col(dflt=0, pos=2),
# 	'time' : Float64Col(dflt=0, pos=3),
# 	'eed' : Float64Col(dflt=0.0, pos=4)
# }

RamaTable = {
	'temp' : Int64Col(dflt=0, pos=0),
	'replicanum' : Int64Col(dflt=0, pos=1),
	'seqnum' : Int64Col(dflt=0, pos=2),
	#'time' : Float64Col(dflt=0, pos=3),
	'phi' : Float64Col(dflt=0.0, pos=4),
	'psi' : Float64Col(dflt=0.0, pos=5),
	'residue' : StringCol(16,pos=6)
}

QTable = {
	'temp' : Int64Col(dflt=0, pos=0),
	'replicanum' :  Int64Col(dflt=0, pos=1),
	'seqnum' : Int64Col(dflt=0, pos=2),
	'framenative' : Int64Col(dflt=0, pos=3),
	'total' : Int64Col(dflt=0, pos=4),
	'nativetotal' : Int64Col(dflt=0, pos=5),
	'Q' : Float64Col(dflt=0, pos=6)
}

ContactMapTable = {
	'temp' : Int64Col(dflt=0, pos=0),
	'replicanum' : Int64Col(dflt=0,pos=1),
	'seqnum' : Int64Col(dflt=0, pos=2),
	'res1' : Int64Col(dflt=0, pos=3), 'res2' : Int64Col(dflt=0, pos=4), 'res3' : Int64Col(dflt=0, pos=5),
	'res4' : Int64Col(dflt=0, pos=6), 'res5' : Int64Col(dflt=0, pos=7), 'res6' : Int64Col(dflt=0, pos=8),
	'res7' : Int64Col(dflt=0, pos=9), 'res8' : Int64Col(dflt=0, pos=10), 'res9' : Int64Col(dflt=0, pos=11),
	'res10' : Int64Col(dflt=0, pos=12), 'res11' : Int64Col(dflt=0, pos=13), 'res12' : Int64Col(dflt=0, pos=14),
	'res13' : Int64Col(dflt=0, pos=15), 'res14' : Int64Col(dflt=0, pos=16), 'res15' : Int64Col(dflt=0, pos=17),
	'res16' : Int64Col(dflt=0, pos=18), 'res17' : Int64Col(dflt=0, pos=19), 'res18' : Int64Col(dflt=0, pos=20),
	'res19' : Int64Col(dflt=0, pos=21), 'res20' : Int64Col(dflt=0, pos=22), 'res21' : Int64Col(dflt=0, pos=23), 
	'res22' : Int64Col(dflt=0, pos=24), 'res23' : Int64Col(dflt=0, pos=25), 'res24' : Int64Col(dflt=0, pos=26), 
	'res25' : Int64Col(dflt=0, pos=27), 'res26' : Int64Col(dflt=0, pos=28), 'res27' : Int64Col(dflt=0, pos=29), 
	'res28' : Int64Col(dflt=0, pos=30), 'res29' : Int64Col(dflt=0, pos=31), 'res30' : Int64Col(dflt=0, pos=64), 
	'res31' : Int64Col(dflt=0, pos=33), 'res64' : Int64Col(dflt=0, pos=34), 'res33' : Int64Col(dflt=0, pos=35), 
	'res34' : Int64Col(dflt=0, pos=36), 'res35' : Int64Col(dflt=0, pos=37), 'res36' : Int64Col(dflt=0, pos=38), 
	'res37' : Int64Col(dflt=0, pos=39), 'res38' : Int64Col(dflt=0, pos=40), 'res39' : Int64Col(dflt=0, pos=41), 
	'res40' : Int64Col(dflt=0, pos=42), 'res41' : Int64Col(dflt=0, pos=43), 'res42' : Int64Col(dflt=0, pos=44), 
	'res43' : Int64Col(dflt=0, pos=45), 'res44' : Int64Col(dflt=0, pos=46), 'res45' : Int64Col(dflt=0, pos=47), 
	'res46' : Int64Col(dflt=0, pos=48), 'res47' : Int64Col(dflt=0, pos=49), 'res48' : Int64Col(dflt=0, pos=50), 
	'res49' : Int64Col(dflt=0, pos=51), 'res50' : Int64Col(dflt=0, pos=51), 'res51' : Int64Col(dflt=0, pos=52), 
	'res52' : Int64Col(dflt=0, pos=53), 'res53' : Int64Col(dflt=0, pos=54), 'res54' : Int64Col(dflt=0, pos=55), 
	'res55' : Int64Col(dflt=0, pos=56), 'res56' : Int64Col(dflt=0, pos=57), 'res57' : Int64Col(dflt=0, pos=58), 
	'res58' : Int64Col(dflt=0, pos=60), 'res59' : Int64Col(dflt=0, pos=61)
}

 
DefaultTable = {
	'time' : Int64Col(pos=0),
	'replicanum' : Int64Col(pos=1),
	'seqnum' : Int64Col(pos=2),
	'temp' : Int64Col(pos=3),
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

