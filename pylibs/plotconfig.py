# -*- coding: utf-8 -*-
import matplotlib.pylab as plt
import ConfigParser

def parseConfig():
	#look for a configuration file in current directory
	#if not found then looking of ~/.analysisrc
	#if not found then quit with error
	#if no config in current dir then plot with default params
	config = ConfigParser.ConfigParser()
	default = '~/config.ini'
	current = os.path.join(os.getcwd(), 'config.ini')
	config.readfp(open(os.path.expanduser(default)))
	config.read(current)

	if not os.path.exists(current):
		print "you did not have a ini file in the current directory,",
	default, "will be used to set the entire analysis settings"

	return config
	
def setConfig(config):
	global analysisParams
	
	if config.get('outputting', 'txt'):
		txtformat = config.get('formatting', 'text')
		analysisParams['txtext'] = txtformat
		
	if config.get('outputting', 'plot'):
		#set plot settings
		if config.get('plotting', 'axes_tight'):
			plt.axis('tight')

		plt.grid(b=config.get('plotting', 'grid'), which='both')	
		plt.rcParams['savefig.dpi'] = config.get('plotting',
'resolution')
		plt.rcParams['font.size'] = config.get('plotting', 'font_size')
		figformat = config.get('formatting', 'figure')
		analysisParams['figext'] = figformat
