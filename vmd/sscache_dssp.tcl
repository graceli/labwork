
proc load_dssp_xpm {} {
	global dssp_xpm
	set infile [open "/home/grace/work/scripts/tcl_test/test.parse" r]
	set frameNum 0
	while { [gets $infile line] >= 0 } {
        	#convert lines read into lists
	        set lst [split $line " "]
        	#store the list in a hash (emulation of a 2D array)
	        set dssp_xpm($frameNum) $lst
        	#puts $dssp_xpm($frameNum)
	        incr frameNum
	}
}

# start the cache for a given molecule
proc start_sscache {{molid top}} {
    load_dssp_xpm

    #global sscache_data
    global dssp_xpm

    if {! [string compare $molid top]} {
	set molid [molinfo top]
    }

    global vmd_frame
    # set a trace to detect when an animation frame changes
    trace variable vmd_frame($molid) w sscache
    return
}

# remove the trace (need one stop for every start)
proc stop_sscache {{molid top}} {
    if {! [string compare $molid top]} {
	set molid [molinfo top]
    }
    global vmd_frame
    trace vdelete vmd_frame($molid) w sscache
    return
}


# reset the whole secondary structure data cache
proc reset_sscache {} {
    global sscache_data
    if [info exists sscache_data] {
        unset sscache_data
    }
    return
}

# when the frame changes, trace calls this function
proc sscache {name index op} {
    # name == vmd_frame
    # index == molecule id of the newly changed frame
    # op == w
    #global sscache_data
    global dssp_xpm
    set sel [atomselect $index "protein name CA"]
    set frame [molinfo $index get frame]
    $sel set structure $dssp_xpm($frame)
    return
}
