#!/bin/sh

# Variables set by the configuration script:
LIBS         = -lxml2 -lmd -lgmx -lnsl -lfftw3f
LDFLAGS      = -L /project/pomes/cneale/GPC/exe/intel/gromacs-4.0.5/exec/lib -L /project/pomes/cneale/GPC/exe/intel/fftw-3.1.2/exec/lib  
CFLAGS       = -O3 -w -funroll-all-loops -std=gnu99  -I/project/pomes/cneale/GPC/exe/intel/gromacs-4.0.5/exec/include -I/project/pomes/cneale/GPC/exe/intel/gromacs-4.0.5/exec/include/gromacs
CC           = cc
LD           = $(CC)

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

all: g_pp_nonpolar g_inositol_residue_nonpolar_v2

g_pp_nonpolar:  g_pp_nonpolar.o
	$(LD) $(LDFLAGS) -o /home/grace/bin/$@ g_pp_nonpolar.o $(LIBS)

g_inositol_residue_nonpolar_v2: g_inositol_residue_nonpolar_v2.o
	$(LD) $(LDFLAGS) -o /home/grace/bin/$@ g_inositol_residue_nonpolar_v2.o $(LIBS)

g_eed: g_eed.o
	$(LD) $(LDFLAGS) -o /home/grace/bin/$@ g_eed.o $(LIBS)

g_mdmat_g: g_mdmat_g.o
	$(LD) $(LDFLAGS) -o /home/grace/bin/$@ g_mdmat_g.o $(LIBS)

clean:
	rm *.o
