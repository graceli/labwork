
#replica object to encapsulate a single xtc analysis

	
def main():
	"""docstring for main"""
	
	#path on disk
	#use /dev/shm
	
	
	#initialize loader object
	aloader = Loader()
	
	# detect all tar files on file system
	# If no STDR_running.tar files found, quit
	
	#copy tar file into /dev/shm
	#attempt to inflate tar in /dev/shm
	
	# analyze xtc files
	for xtcfile in xtcList:
		#get basename
		#pass in path/to/file/ and basename
		xtc = Xtc(path, xtcfile, tprfile)
		
		rgpath = xtc.rg()
		aloader.load('rg', rowtypes.RGTable, 4, rgpath)
	
		sasapath = xtc.sasa()
		aloader.load('sas', rowtypes.SASTable, 4, sasapath)
		
		#aloader.load('sas', rowtypes.SASTable, replicaMeta)
		#aloader.load('eed', rowtypes.EEDTable, replicaMeta)
		#aloader.load('dihedral', rowtypes.DihedralTable, replicaMeta)
		#aloader.load('energy', rowtypes.EnergyTable, replicaMeta)
		
	

if __name__ == '__main__':
	main()