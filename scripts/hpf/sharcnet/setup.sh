#make sure to have an existing gro with molname

if [ -z $1 ]; then
	echo "usage: setup.sh <molname>";
	exit;
fi
molname=$1;

#makes an ndx file
echo "q" | make_ndx -f ${molname}.gro -o ${molname}.ndx

#makes an md0 dir and put the initial gro file into the md0_success folder
mkdir md0_success
cp ${molname}.gro md0_success/${molname}_md0_deshuffleddesorted.gro

#echo a reminder to have an mdp file and modify head.sh appropriately
echo 
echo 
echo "####### IMPORTANT: Please also place ${molname}.mdp and modify head.sh appropriately #######"

