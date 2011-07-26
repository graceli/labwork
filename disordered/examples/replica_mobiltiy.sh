# Hi Grace,
# 
# Actually, replica mobility reports on the decorrelation time of the potential energy, but a hidden barrier is just a special case of this. A local minimum on the output from this script indicates a region of longer than average decorrelation times for the potential energy.
# 
# here is the script to analyze replica mobility. I have copied it to this email, but you should probably copy it from scinet:/scratch/cneale/GPC/ARG/DR_as_indentification_tool/forgrace.sh to ensure that the email doesn't split lines or something.
# 
# Note that you need to change the id to your chosen 2-letter code. Also, the script expects a single forcedatabase and I think that you have many. I suggest starting by just using one forcedatabase that is fairly long. Later, you can create many database.txt files and put them together in the proper order to run the entire thing through this script.
# 
# Good luck, and I'd love to see the resulting plot once you have run it.
# 
# Sarah and David, if you are interested, we can talk about how to do this with your data, which I assume is in a different format and can not simply be analyzed by these programs. Still, it is a very simple analysis.
# 
# Thanks,
# Chris.

id=ta

DOSTART=1
DOMORE=1

SHIFTA=3
SHIFTB=5
SHIFTC=8

if((DOSTART)); then

mkdir -p UP
cd UP
if((DOMORE)); then
 cp ../${id}.forcedatabase .
 /project/pomes/cneale/GPC/exe/DR/DR_version2.3.2/bin/analyse_force_database ${id}.script -d database.txt -t 0 > analysis.ps
fi
/project/pomes/cneale/GPC/exe/DR/DR_version2.4.1/source/extractDatabase_mobility3 database.txt ${SHIFTA}>za
/project/pomes/cneale/GPC/exe/DR/DR_version2.4.1/source/extractDatabase_mobility3 database.txt ${SHIFTB}>zb
/project/pomes/cneale/GPC/exe/DR/DR_version2.4.1/source/extractDatabase_mobility3 database.txt ${SHIFTC}>zc
cd ../

module load graphics
gnuplot  << EOF
set term post
set out 'a.ps'
plot 'za' u 1:2 w boxes lt 1, 'zb' u 1:2 w boxes lt 1 fs solid 0.3, 'zc' u 1:2 w boxes lt 1 fs solid 0.6
EOF