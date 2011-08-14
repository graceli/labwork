proc draw_plane {molid color} {

	set x [molinfo $molid get a]
	set y [molinfo $molid get b]
	set z [molinfo $molid get c]

	
	draw color $color 
	draw line "0 0 $z/2" "$x 0 $z/2" 
	draw line "0 0 $z/2" "0 $y $z/2"
	draw line "$x 0 $z/2" "$x $y $z/2"
	draw line "0 $y $z/2" "$x $y $z/2"
}

