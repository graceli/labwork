proc draw_box {molid color} {

	set x [molinfo $molid get a]
	set y [molinfo $molid get b]
	set z [molinfo $molid get c]
	
	draw color $color 
	draw line "0 0 0" "$x 0 0"
	draw line "0 0 0" "0 $y 0"
	draw line "0 0 0" "0 0 $z"
	draw line "$x 0 0" "$x $y 0"
	draw line "$x 0 0" "$x 0 $z"
	draw line "0 $y 0" "$x $y 0"
	draw line "0 $y 0" "0 $y $z"
	draw line "0 0 $z" "0 $y $z"
	draw line "0 0 $z" "$x 0 $z"
	draw line "$x $y 0" "$x $y $z"
	draw line "0 $y $z" "$x $y $z"
	draw line "$x 0 $z" "$x $y $z"
}

