#!/bin/sh

for file in `ls *water*sc.xvg`; do 
	echo "replacing $file"
	sed --in-place=.bck -e 's/@/#/g' $file 
done
