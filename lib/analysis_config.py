import pylab

STARTFRAME = 20000
LASTFRAME = 90000
NMOLECULES = 4
RATIO = { "15to4" : "4:1", "45to4" : "10:1" }
HEADER_NONPOLAR = "# intermolecular nonpolar contacts per peptide"
RUNTIME_NS = 180
LINE_COLOR = { 'scyllo' : 'red', 'chiro' : 'blue', 'water' : 'green', 'control': '#000000' }
SHADED_COLOR = { 'scyllo' : 'lightblue', 'chiro' : 'pink', 'water' : 'chartreuse' }
LABEL = { 'scyllo' : 'scyllo', 'chiro' : 'chiro', 'water' : 'no inositol' }
ISOMER_LIST = ['scyllo', 'chiro', 'water']


def plot_settings(setting="default"):
	pylab.rc('text', usetex=True)
	if setting == "pub":
		pylab.rcParams['figure.figsize'] = [3,3]
		pylab.rcParams['font.size'] = 9
		pylab.rcParams['legend.fontsize'] = 'small'
	else:
		pylab.rcParams['figure.figsize'] = [8,6]
		pylab.rcParams['font.size'] = 12