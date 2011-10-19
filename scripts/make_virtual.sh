#!/bin/sh

#set up python virtual environment
python python/virtualenv-1.5.1/virtualenv ENV
# note to change this /home/grace/.pydistutils.cfg
# source ~/ENV/bin/activate
# deactivate 

#install GromacsWrapper and test it
easy_install -f http://sbcb.bioch.ox.ac.uk/oliver/download/Python GromacsWrapper
python -c "import gromacs"

#note that ipython doesn't work for some reason with easy_install under virtualenv
