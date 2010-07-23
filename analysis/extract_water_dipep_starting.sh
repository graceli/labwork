#!/bin/bash

echo "0" | trjconv -f inos_ala_dipep_md_nvt.xtc -o water_dipep_starts/water_dipep4.pdb -s inos_ala_dipep_md_nvt.tpr -b 4 -e 5
# this is deprecated
#echo "0" | trjconv -f inos_ala_dipep_md_nvt.xtc -o water_dipep_starts/water_dipep7.pdb -s inos_ala_dipep_md_nvt.tpr -b 7 -e 8
# and replaced with this:
echo "0" | trjconv -f inos_ala_dipep_md_nvt.xtc -o water_dipep_starts/water_dipep49.pdb -s inos_ala_dipep_md_nvt.tpr -b 49 -e 50 
echo "0" | trjconv -f inos_ala_dipep_md_nvt.xtc -o water_dipep_starts/water_dipep93.pdb -s inos_ala_dipep_md_nvt.tpr -b 93 -e 94
echo "0" | trjconv -f inos_ala_dipep_md_nvt.xtc -o water_dipep_starts/water_dipep108.pdb -s inos_ala_dipep_md_nvt.tpr -b 108 -e 109
echo "0" | trjconv -f inos_ala_dipep_md_nvt.xtc -o water_dipep_starts/water_dipep189.pdb -s inos_ala_dipep_md_nvt.tpr -b 189 -e 190
echo "0" | trjconv -f inos_ala_dipep_md_nvt.xtc -o water_dipep_starts/water_dipep203.pdb -s inos_ala_dipep_md_nvt.tpr -b 203 -e 204
