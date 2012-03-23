#!/bin/sh
#PBS -q qwork@mp2 
#PBS -l walltime=00:24:00 -l nodes=2:ppn=1
#PBS -N analysis

# Extract protein, inositol, and solvent Z density profiles
# Want to see how each of these components are distributed around the system wrt the bilayer.
# Anything of these things partitioning to or into the bilayer?  hard to imagine but not impossible.

# TODO add index groups
set -x
set -e 
set -u

#TEST="-dt 1000"
TEST=""

trap "clean; exit $?" TERM INT SIGINT EXIT

#cd $PBS_O_WORKDIR

function clean {
	mkdir -p analysis/density
	cp -rp /dev/shm/grace/* analysis/density/
	rm -rf /dev/shm/grace
}

function plot {
for ID in `seq 0 43`; do
cd $ID
gnuplot <<EOF
set term png
set output "/dev/shm/grace/densities_${ID}.png"

set datafile commentschars "#@"
plot "protein_density_${i}.xvg" w l, "sol_density_${i}.xvg" w l, "inositol_density_${i}.xvg" w l, "bilayer_density_${i}.xvg" w l, "headgroup_density_${i}.xvg" w l
EOF
done
cd ../
}

function density {
	INDEX_GROUP=$1	
	for i in `seq 0 43`; do
		cd $i
	if [ ! -e "toxicity.ndx" ]; then
make_ndx_mpi -f prod.tpr -o toxicity.ndx <<EOF
a P*
name 18 headgroups
q
EOF
	fi

		echo INDEX_GROUP | g_density_mpi -f prod.part0003.xtc -s prod.tpr -o /dev/shm/grace/${INDEX_GROUP}_density_${i}.xvg -n toxicity.ndx $TEST 2> /dev/shm/grace/${INDEX_GROUP}_density_${i}.out >&2 &
		cd ../
	done
}

wait

density "Protein"
wait
density "SOL"
wait
density "INS"
wait
density "headgroups"
wait
density "POP"
wait

plot

