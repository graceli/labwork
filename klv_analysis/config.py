import pylab
import plot_settings as ps

LASTFRAME = 90000
NMOLECULES = 4
RATIO = { "15to4" : "4:1", "45to4" : "10:1" }
HEADER_NONPOLAR = "# intermolecular nonpolar contacts per peptide"
RUNTIME_NS = 180
LINE_COLOR = { 'scyllo' : 'blue', 'chiro' : 'red', 'water' : 'green' }
SHADED_COLOR = { 'scyllo' : 'lightblue', 'chiro' : 'pink', 'water' : 'chartreuse' }
LABEL = { 'scyllo' : 'scyllo', 'chiro' : 'chiro', 'water' : 'no inositol' }

# plot_settings thirdparty
# config.configure_plot()
# config.ps.set_mode("publish")

def configure_plot():
	print "set plot configuration"
	# pylab.rcParams['backend'] = 
	# pylab.clf()
	pylab.rcParams['text.usetex'] = False
	pylab.rcParams['font.size'] = 10
	pylab.rcParams['font.family'] = 'Arial'
	# pylab.rcParams['legend.fontsize']='small'
	pylab.rcParams['legend.fontsize'] = 8
	#subplots_adjust(left=0, bottom=0, right=0.001, top=0.001, wspace=0, hspace=0)
	pylab.rcParams['figure.figsize'] = [3.5,3.5]

