module load graphics
for i in 1 2 3; do
  /project/pomes/cneale/GPC/exe/DR/DR_version2.1.8_extended/bin/analyse_force_database_nofirstpage 00.script -a $i > analysis${i}.ps
  ps2pdf analysis${i}.ps analysis${i}.pdf
done
