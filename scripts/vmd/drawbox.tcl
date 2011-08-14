# draws the box with dimensions as specified by the structure file
# usage example: "draw_box top white" draws boxes with white lines around 
# your system
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
	
	puts "Box dimensions are $x x $y x $z A^3"
}

# same as drawbox but allows to specify an arbitrary box dimension( in nanometers)
# eg: draw_box_size top 3 3 3 white
# draws a box with length 3 nm on each side
proc draw_box_size {molid x y z color} {
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

#draws a box with an arbitrary origin
proc draw_box_origin {molid o1 o2 o3 x y z color} {
	draw color $color
	set xs [expr $o1+$x];
	set ys [expr $o2+$y];
	set zs [expr $o3+$z];

	draw line "$o1 $o2 $o3" "$xs $o2 $o3"
	draw line "$o1 $o2 $o3" "$o1 $ys $o3"
	draw line "$o1 $ys $o3" "$xs $ys $o3"
	draw line "$xs $ys $o3" "$xs $o2 $o3"

	draw line "$o1 $o2 $zs" "$xs $o2 $zs"
	draw line "$o1 $o2 $zs" "$o1 $ys $zs"
	draw line "$o1 $ys $zs" "$xs $ys $zs"
	draw line "$xs $ys $zs" "$xs $o2 $zs"

	draw line "$o1 $ys $o3" "$o1 $ys $zs"
	draw line "$xs $ys $o3" "$xs $ys $zs"
	draw line "$o1 $o2 $o3" "$o1 $o2 $zs"
	draw line "$xs $o2 $o3" "$xs $o2 $zs"
}
