# generic 2d plotting script for gnuplot
set term post color enhanced solid 8 portrait
set output "ala_scy_bidentate_overlay.ps"

unset surface
set view 0,0
#set hidden3d
set zrange [0.01:10]
set contour base 
set cntrparam linear
set cntrparam levels incremental 0.001,0.3,3.5

set tics in
#set tics level 0.5
#set tic scale 1 0.5
set grid
set key right box title "kcal/mol" 

set size square
set size 1,1
set origin 0,0
set xlabel "phi"
set ylabel "psi"

#compress contour plotted in 3D format into a 2D format via table terminal
splot [-180:180][-180:180] 'awat_5deg.pmf' w l
set table 'awat_contour.dat'
replot
unset table

set multiplot layout 4,3 rowsfirst title "alanine dipeptide scyllo inositol PMFs"
unset key
#draw contour in 2D with the (phi,psi) rama map
plot [-180:180][-180:180] 'ala_scy_CO0CO2_12.dat' , 'awat_contour.dat' w l
plot [-180:180][-180:180] 'ala_scy_CO0CO2_13.dat' , 'awat_contour.dat' w l
unset grid
plot [0:0.1][0:0.1] 0 w dots
set grid
plot [-180:180][-180:180] 'ala_scy_NH1CO2_12.dat' , 'awat_contour.dat' w l
unset grid
plot [0:0.1][0:0.1] 0 w dots
set grid
unset grid
plot [0:0.1][0:0.1] 0 w dots
set grid
plot [-180:180][-180:180] 'ala_scy_CO0NH3_12.dat' , 'awat_contour.dat' w l
plot [-180:180][-180:180] 'ala_scy_CO0NH3_13.dat' , 'awat_contour.dat' w l
unset grid
plot [0:0.1][0:0.1] 0 w dots
set grid
plot [-180:180][-180:180] 'ala_scy_NH1NH3_12.dat' , 'awat_contour.dat' w l
unset grid
plot [0:0.1][0:0.1] 0 w dots
set grid
unset grid
plot [0:0.1][0:0.1] 0 w dots
set grid


unset multiplot
