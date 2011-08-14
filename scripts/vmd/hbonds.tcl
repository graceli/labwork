proc hb {} {
	mol representation "HBonds 3.5 60 5"
	mol selection "resname INS or backbone and not name \"CA*\""
	mol addrep top
}

