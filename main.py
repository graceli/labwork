# -*- coding: utf-8 -*-
from loader import *

agent = Loader()

agent.load('rg', rowtypes.RGTable,'data/rg')
agent.load('sas',rowtypes.SASTable, 'data/sasa')

