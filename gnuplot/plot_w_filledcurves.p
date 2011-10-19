set term post color enhanced solid
set output "klv_monomer_eed_w_chiro_scyllo_noins.eps"

set xlabel "End-to-end distance (A)"
set ylabel "Probability"

#plot [4:21] 'scyllo_eed_withError.histo' w l lw 2 t 'scyllo-', '' u 1:2:3 w yerrorbars t '', 'chiro_withError_0.2ang.histo1000' w l lw 2 t 'chiro-', '' u 1:2:3 w yerrorbars t '', 'water_eed_withError_0.2ang.histo1000' w l lw 2 t 'no inositol', '' u 1:2:3 w yerrorbars t ''

plot [4:21] 'water_eed_withError_0.2ang.histo1000'  u 1:($2-$3):($2+$3) w filledcurves lc rgb "#dddddd" t '', 'scyllo_eed_withError.histo' u 1:($2-$3):($2+$3) w filledcurves lc rgb "#dddddd" t '', 'chiro_withError_0.2ang.histo1000'  u 1:($2-$3):($2+$3) w filledcurves lc rgb "#dddddd" t '', 'water_eed_withError_0.2ang.histo1000'  u 1:2 w l t 'no inositol', 'scyllo_eed_withError.histo' u 1:2 w l t 'scyllo-inositol', 'chiro_withError_0.2ang.histo1000'  u 1:2 w l t 'chiro-inositol'
