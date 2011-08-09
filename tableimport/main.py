#!/usr/bin/env python
from loader import *
import os

nfixedcols = 4
agent = Loader()

rgpath = os.path.join(os.getcwd(), 'data/rg')
sasapath = os.path.join(os.getcwd(), 'data/sasa')

agent.load('rg', rowtypes.RGTable, nfixedcols, rgpath)
agent.load('sas',rowtypes.SASTable, nfixedcols,  sasapath)

