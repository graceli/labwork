#!/bin/bash

#PBS -l walltime=48:00:00
#PBS -q archive 
#PBS -N hpss_index
#PBS -j oe
#PBS -m e

/scinet/gpc/bin/ish hindex

