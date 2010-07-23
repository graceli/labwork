proc draw_ins {} {
mol representation "hbonds"
mol selection "backbone and not name \"CA*\""
mol addrep top
mol representation "vdw"
mol selection "resname INS"
mol addrep top
mol representation "NewCartoon"
mol selection "all and not resid 2 to 9"
mol addrep top
mol representation "NewCartoon"
mol selection "resid 2 to 9"
mol color "colorid 4"
mol addrep top
}

