#!/bin/sh

echo "This script makes a nosol.gro and a nosol.ndx file";

if [ ! -e *.gro ]; then
        echo "ERROR: gro file does not exist";
else
        grofile=`ls *.gro`;
        echo -e "!\"SOL\"\nq" | make_ndx -f $grofile -o nosol.ndx
	echo "nosol.ndx is successfully made";

        echo -e "!SOL" | trjconv -f $grofile -s $grofile -n nosol.ndx -o nosol.gro
	echo "nosol.gro is successfully outputted";
fi

