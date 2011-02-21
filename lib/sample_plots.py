#!/usr/bin/python

# Author: Varun Hiremath <vh63@cornell.edu>
# Created: Thu,  2 Apr 2009 05:21:30 -0400

import scipy, pylab
from scipy import io
import plot_settings as ps

error = io.read_array("error.op")
nrs = error[:,0]

# Publishable quality image
ps.set_mode("publish")
pylab.figure(1)
pylab.axes([0.125,0.2,0.95-0.125,0.95-0.3])

for i in range(6, len(error[0])):
    pylab.plot(nrs, error[:,i], ps.lps[i], label="$\Phi_{"+str(i-5)+"}$")

pylab.legend(numpoints=1)
pylab.xlabel("$n_{rs}$")
pylab.ylabel("error ($\epsilon$)")
pylab.title("Plot of error")
pylab.savefig("error_publish.eps")
pylab.savefig("error_publish.png")

# Medium size image
ps.set_mode("medium")
pylab.figure(2)
# pylab.axes([0.125,0.2,0.95-0.125,0.95-0.3])

for i in range(6, len(error[0])):
    pylab.plot(nrs, error[:,i], ps.lps[i], label="$\Phi_{"+str(i-5)+"}$")

pylab.legend(numpoints=1)
pylab.xlabel("$n_{rs}$")
pylab.ylabel("error ($\epsilon$)")
pylab.title("Plot of error")
pylab.savefig("error_medium.eps")
pylab.savefig("error_medium.png")
