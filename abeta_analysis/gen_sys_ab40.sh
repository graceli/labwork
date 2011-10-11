#!/bin/sh

function clean {
	rm \#*
	rm -rf water scyllo_* chiro_*
}

function error_exit {
	if [ "$?" > "0" ]; then
		echo "Error!"
		exit 1
	fi
} 

clean;

CONF="conf2.gro"
TOP="abeta40.top"
TOP_INS="start/abeta40_inositol.top"
start=1
end=100
# water

mkdir water
cd water
for i in `seq $start $end`; do 
	mkdir sys${i}
	cp ../start/$TOP sys${i}/abeta40_water.top
	top="abeta40_water.top"
	cp ../start/$CONF sys${i}
	cp ../start/*.mdp sys${i}
	cp ../start/*.itp sys${i}
	cp ../start/*.dat sys${i}
	sed -e "s/INDEX/${i}/g" -e "s/TOP/$top/g" ../start/sys0.sh > sys$i/sys${i}.sh
	cd sys${i}
	genbox -cp $CONF -cs ~/gro/tip3.gro -p abeta40_water.top -o abeta40_water.gro
	touch empty.mdp
	grompp -f empty.mdp -c abeta40_water.gro -p abeta40_water.top -o out1.tpr
	echo 12 | genion -s out1.tpr -conc 0.01 -nname CL- -pname NA+ -o ions_water.gro -p abeta40_water.top
	grompp -f empty.mdp -c ions_water.gro -p abeta40_water.top -o ions.tpr
	echo 12 | genion -s ions.tpr -nn 12 -np 0 -nname CL- -pname NA+ -o ions_counter_water.gro -p abeta40_water.top
	
	#grompp
	grompp -f em.mdp -c ions_counter_water.gro -p abeta40_water.top -o em.tpr
	cd ../
done
cd ../

for iso in scyllo chiro; do 
    for ratio in 64 128; do
		mkdir ${iso}_${ratio}
		cd ${iso}_${ratio}
		for i in `seq $start $end`; do 	
			mkdir sys${i}
			cp ../start/abeta40_inositol_${ratio}.top sys${i}/abeta40_inositol_${ratio}.top
			top="abeta40_inositol_${ratio}.top"
			cp ../start/conf2.gro sys${i}
			cp ../start/*.mdp sys${i}
			cp ../start/*.itp sys${i}
			cp ../start/*.dat sys${i}
			sed -e "s/INDEX/${i}/g" -e "s/TOP/$top/g" ../start/sys0.sh > sys${i}/sys${i}.sh
			cd sys${i}
			genbox -cp $CONF -ci ~/gro/${iso}_em.gro -nmol $ratio -o abeta40_${iso}_${ratio}_box.gro -p abeta40_inositol_${ratio}.top
			genbox -cp abeta40_${iso}_${ratio}_box.gro -cs ~/gro/tip3.gro -p abeta40_inositol_${ratio}.top -o abeta40_${iso}_${ratio}.gro
			touch empty.mdp
			grompp -f empty.mdp -c abeta40_${iso}_${ratio}.gro -p abeta40_inositol_${ratio}.top -o out1.tpr
			echo 13 | genion -s out1.tpr -conc 0.01 -nname CL- -pname NA+ -o ions.gro -p abeta40_inositol_${ratio}.top
			grompp -f empty.mdp -c ions.gro -p abeta40_inositol_${ratio}.top -o ions.tpr
			echo 13 | genion -s ions.tpr -nn 12 -np 0 -nname CL- -pname NA+ -o ions_counter.gro -p abeta40_inositol_${ratio}.top
			# grompp
			grompp -f em.mdp -c ions_counter.gro -p abeta40_inositol_${ratio}.top -o em.tpr
			cd ../
		done
		cd ../
    done
done


