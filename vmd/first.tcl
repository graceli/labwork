set file "test.parse"
set infile [open $file r]
set frameNum 0
while { [gets $infile line] >= 0 } {
	#convert lines read into lists
	set lst [split $line " "]
	#store the list in a hash (emulation of a 2D array)
	set name($frameNum) $lst
	puts $name($frameNum)
	
	incr frameNum
}
