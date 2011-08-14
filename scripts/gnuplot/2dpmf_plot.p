# generic 2d plotting script for gnuplot

# for setting the color map view
#set view map
#set style data pm3d

# do no show surface
unset surface

#sets contours and parameters
set view 1,0,1,
set zrange [-8:1]
set contour base
set cntrparam linear
set cntrparam levels incremental -6, 0.3, -2

set tics in 
set ticslevel 0.5
set ticscale 1 0.5
set grid

set origin 0.0,0.0
splot [-180:180][-180:180] 'plot' w l
