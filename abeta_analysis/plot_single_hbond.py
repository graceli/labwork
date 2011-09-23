from pylab import *
import numpy
import utils

#def smooth(data, window_size, time_present=True, timestep=1):
rcParams['xtick.labelsize']='10'
rcParams['ytick.labelsize']='12'
rcParams['legend.fontsize']='12'
rcParams['figure.figsize'] = [2.5,2.5]
rcParams["axes.titlesize"]='small'

A=utils.smooth(numpy.genfromtxt('chain_0_1_hbonds.xvg', skip_header=20), 1000, timestep=2)
B=utils.smooth(numpy.genfromtxt('chain_1_2_hbonds.xvg', skip_header=20), 1000, timestep=2) 
C=utils.smooth(numpy.genfromtxt('chain_2_3_hbonds.xvg', skip_header=20), 1000, timestep=2)
D=utils.smooth(numpy.genfromtxt('chain_3_4_hbonds.xvg', skip_header=20), 1000, timestep=2)

print A

time = A[:,0]/1000.0
plot(time, A[:,1], label="chain 1-2")
plot(time, B[:,1], label="chain 2-3")
plot(time, C[:,1], label="chain 3-4")
plot(time, D[:,1], label="chain 4-5")
ylim(0, 30)
xlim(0, 140)
grid(True)
# xlabel('Time (ns)')
# ylabel('Number of interchain hydrogen bonds')
#legend(loc='lower right')
savefig('chain_hbonds.png')


