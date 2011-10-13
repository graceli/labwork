import glob

def main():
	isomer = ['glycerol', 'chiro', 'scyllo', 'water']
	ratios = [ 15, 64 ]
	for iso in isomer:
		for ratio in ratios:
			xtc_list = glob.glob("*%(iso)s*%(ratio)d*.xtc" % vars())
			for file in sorted(xtc_list):
				print file
			print "total:", len(xtc_list)



if	__name__=='__main__':
	main()
