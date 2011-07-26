#!/usr/bin/python

# generate a list of tar files to process
# persist this in some sort of fast python readable file (pickle?)
# file => processed?
# lgwXXX => true
# lgwYYY => false

#set -x
set -e
set -u
set -o noclobber

DATA=/scratch/grace/drSH3/PRIOR_TO_RESTART_Fri_Dec_24_17:32:02_EST_2010/output/data
SHM=/dev/shm

for file in `ls $DATA/*tmp*`; do
    echo "processing $file ..."
    tar xf $file -C $SHM
    cd $SHM
    if [ ! -e 300 ]; then
        mkdir 300
    fi

    for xtc in `ls *.xtc`; do
        echo $xtc
        bname=`basename $xtc .xtc`
        edrfile="${bname}.edr"
        echo Temperature | g_energy -f $edrfile
        T=`echo Temperature | g_energy -f $edrfile 2> /dev/null | grep Temperature | awk '{print $2}'`
        echo $bname
        echo $edrfile
        echo $T
        exit 0
    done
done