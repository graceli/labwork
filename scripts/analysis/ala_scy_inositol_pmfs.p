# generic 2d plotting script for gnuplot
set term post color enhanced solid 8 portrait
set output "ala_scy_pmfs_0.1.ps"

unset surface
set view 0,0
set hidden3d
set zrange [0.01:10]
set contour base 
set cntrparam linear
set cntrparam levels incremental 0.001,0.1,0.8

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
set multiplot layout 3,1 title "alanine dipeptide scyllo inositol PMFs"
set title "ala-scyllo non-binding"
set rmargin 0
splot [-180:180][-180:180] 'ala_scy_none_5deg.pmf' w l  
set title "ala-scyllo mono"
splot [-180:180][-180:180] 'ala_scy_mono_5deg.pmf' w l  
set title "ala-scyllo bidentate"
splot [-180:180][-180:180] 'ala_scy_bidentate_5deg.pmf' w l  
unset multiplot
