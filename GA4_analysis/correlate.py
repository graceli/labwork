import tables

def main():
	"""docstring for main"""
	NINS = 4
	
	angles = tables.openFile('GA4_angles_analysis.h5')
	contacts = tables.openFile('GA4_beta_analysis.h5')
	
	f60 = open('angles_vs_contacts_theta60.txt', 'w')
	f30_60 = open('angles_vs_contacts_theta30-60.txt', 'w')
	f30 = open('angles_vs_contacts_theta30.txt', 'w')
	
	for iso in ['chiro', 'scyllo']:
		for i in range(0, 3):
			if i == 1 and iso == 'chiro':
				continue
				
			angles_table = angles.getNode('/angle/inf_perfect_%(iso)s_sys%(i)s' % vars())
			polar_contacts_table = contacts.getNode('/polar-contacts-per-inositol/inf_perfect_%(iso)s_sys%(i)s' % vars())
			nonpolar_contacts_table = contacts.getNode('/nonpolar_per_inositol/inf_perfect_%(iso)s_sys%(i)s' % vars())
			
			for ino in range(1, 5):
				index = 2 + (ino - 1) * 3
				for row in angles_table:
					angle = row[ 'col%(index)d' % vars() ]
					theta = angle
					if angle > 90:
						theta = 180 - angle
				
					row_p = polar_contacts_table[row.nrow]
					row_np = nonpolar_contacts_table[row.nrow]	
	
					if theta > 60:
						print >> f60, row.nrow, iso, "sys%(i)d" % vars(), ino, index, theta, row_p[ 'col%(ino)d' % vars() ], row_np[ 'col%(ino)d' % vars() ]
					elif theta > 30 and theta <= 60:
						print >> f30_60, row.nrow, iso, "sys%(i)d" % vars(), ino, index, theta, row_p[ 'col%(ino)d' % vars() ], row_np[ 'col%(ino)d' % vars() ] 
					else:
						print >> f30, row.nrow, iso, "sys%(i)d" % vars(), ino, index, theta, row_p[ 'col%(ino)d' % vars() ], row_np[ 'col%(ino)d' % vars()]
				
			
if __name__ == '__main__':
	main()