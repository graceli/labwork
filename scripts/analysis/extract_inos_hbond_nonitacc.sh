#!/bin/bash
# $1 = xtc file to read from
# $2 = unique id so that previous hbondnum files aren't overwritten

echo "7 15" | g_hbond -f $1 -s inos_ala_dipep_md_nvt.tpr -n inos_hbond_index.ndx -num inos4_hbnum$2.xvg -nonitacc
echo "7 16" | g_hbond -f $1 -s inos_ala_dipep_md_nvt.tpr -n inos_hbond_index.ndx -num inos5_hbnum$2.xvg -nonitacc
echo "7 17" | g_hbond -f $1 -s inos_ala_dipep_md_nvt.tpr -n inos_hbond_index.ndx -num inos6_hbnum$2.xvg -nonitacc
echo "7 18" | g_hbond -f $1 -s inos_ala_dipep_md_nvt.tpr -n inos_hbond_index.ndx -num inos7_hbnum$2.xvg -nonitacc

