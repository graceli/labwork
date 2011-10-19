#!/usr/env/python
import gromacs
#example code for gromacs wrapper
#run gmxcheck and send stderr to /dev/null
gromacs.gmxcheck(f='GA4_oct_sc0_nosol.xtc', stderr=open(os.devnull,'w'))
