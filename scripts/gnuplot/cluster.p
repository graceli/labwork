#set terminal png transparent nocrop enhanced font arial 8 size 500,350 
set terminal png nocrop enhanced font arial 8 size 500,350 
set output 'histograms.2.png'
set bar 1.000000
set boxwidth 0.9 absolute
#set style fill solid 1.00 border -1
set style fill solid 1.00 noborder

#set style rectangle back fc lt -3 fillstyle  solid 1.00 border -1
set style rectangle back fc lt -3 fillstyle solid 1.00 noborder
set key inside right top vertical Right noreverse enhanced autotitles columnhead nobox
set style histogram clustered gap 2 title  offset character 0, 0, 0
#set datafile missing '-'
set style data histograms

set boxwidth 0.7 relative

plot 'scyllo/one_contact_hist.hist' u 2 ti col t 'scyllo', 'chiro/one_contact_hist.hist' u 2 ti col t 'chiro'
