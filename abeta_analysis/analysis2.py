import glob
import os


def comment_xvgr():
	# comment out the xvgr header lines for the chain chain hbond analysis
	for i in range(1,11):
		for chain in range(0,4):
			next_chain = chain + 1
			os.system("sed -e 's/@/#/g' %(i)d/chain_%(chain)d_%(next_chain)d_hbonds.xvg > chch_hbonds/sys%(i)d_chain_%(chain)d_%(next_chain)d_hbonds.dat" % vars())

	# delete the data files for very short (failed) simulation trajs
	files = glob.glob('chch_hbonds/*.dat')
	for f in files:
		size = os.path.getsize(f)/(1024.0*1024)
		if size < 1:
			os.system('rm %(f)s' % vars())
			
	
if __name__ == '__main__':
	main()