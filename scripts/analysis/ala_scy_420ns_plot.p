set size square
set term post color enhanced 10
set output "ala_scy_420ns_plot.ps"
set xrange [-180:180]
set yrange [-180:180]
set key bottom
set multiplot
set size 0.33,0.5
set origin 0, 0
plot "ala_scy_CO0CO2_12.dat" w dots
set origin 0.33, 0
plot "ala_scy_CO0CO2_13.dat" w p
set origin 0.66, 0
plot "ala_scy_CO0NH3_12.dat" w dots
set origin 0, 0.5
plot "ala_scy_CO0NH3_13.dat" w p
set origin 0.33, 0.5
plot "ala_scy_mono.dat" w dots
set origin 0.66, 0.5
plot "ala_scy_NH1CO2_12.dat" w dots
unset multiplot
set multiplot
set origin 0, 0
plot "ala_scy_NH1NH3_12.dat" w dots
set origin 0.33, 0
plot "ala_scy_none.dat" w dots
unset multiplot
