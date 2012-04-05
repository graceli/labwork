#!/bin/sh
#PBS -q qwork@mp2 
#PBS -l walltime=24:00:00 -l nodes=2:ppn=1
#PBS -N analysis

# Extract protein, inositol, and solvent Z density profiles
# Want to see how each of these components are distributed around the system wrt the bilayer.
# Anything of these things partitioning to or into the bilayer?  hard to imagine but not impossible.

# TODO add index groups
set -x
set -u

#TEST="-dt 1000"
TEST=""

trap "clean; exit $?" TERM INT SIGINT EXIT

cd $PBS_O_WORKDIR

base_dir=$PWD

mkdir /dev/shm/grace

function clean {
	mkdir -p $base_dir/analysis/density
	cp -rp /dev/shm/grace/* $base_dir/analysis/density/
	rm -rf /dev/shm/grace
}

function plot {
cd /dev/shm/grace

for i in `seq 0 43`; do
gnuplot <<EOF
set term png
set output "/dev/shm/grace/densities_${i}.png"

set datafile commentschars "#@"
plot "Protein_density_${i}.xvg" w l, "SOL_density_${i}.xvg" w l, "INS_density_${i}.xvg" w l, "POP_density_${i}.xvg" w l, "headgroups_density_${i}.xvg" w l
EOF
done

cd ../
}

function plot_no_inositol {
cd /dev/shm/grace

for i in `seq 0 43`; do
gnuplot <<EOF
set term png
set output "/dev/shm/grace/densities_${i}.png"

set datafile commentschars "#@"
plot "Protein_density_${i}.xvg" w l, "SOL_density_${i}.xvg" w l, "POP_density_${i}.xvg" w l, "headgroups_density_${i}.xvg" w l
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

		echo "Getting the density of $INDEX_GROUP for system=$i ..."
		echo $INDEX_GROUP | g_density_mpi -f prod.xtc -s prod.tpr -o /dev/shm/grace/${INDEX_GROUP}_density_${i}.xvg -n toxicity.ndx $TEST 2> /dev/shm/grace/${INDEX_GROUP}_density_${i}.out >&2 &
		cd ../
	done
}

wait

density "Protein"
wait
density "SOL"
wait

#	density "INS"

#wait
density "headgroups"
wait
density "POP"
wait

plot_no_inositol

