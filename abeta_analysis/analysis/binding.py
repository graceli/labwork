def intersect(nonpolar, polar):
    # get polar matrix
	nrows, ncols = nonpolar.shape

	print "computing the intersecting"
	print "nonpolar has dimensions ", nonpolar.shape
	print nonpolar

	print "polar has dimensions", polar.shape
	print nonpolar
	print polar
	print "polar ", polar[0], polar[0][0]

	p_and_np_count  = 0
	nonpolar_count = 0
	polar_count = 0
	total_inositol = 0
	for r in range(0, nrows):
		for i in range(0, ncols):
			# print nonpolar_new[r][i], polar[r][i]
			if nonpolar[r][i] > 0.0 and polar[r][i] > 0.0:
				p_and_np_count += 1
			elif nonpolar[r][i] > 0.0 and polar[r][i] == 0.0:
				nonpolar_count += 1
			elif nonpolar[r][i] == 0.0 and polar[r][i] > 0.0:
				polar_count += 1
			total_inositol += 1

	print "system", sys, polar_count, nonpolar_count, p_and_np_count, total_inositol


if __name__ == '__main__':
    for ratio in config.ratio_list:
        for isomer in config.isomer_list:
            nonpolar_h5file = tables.openFile("nonpolar_contacts_" + str(isomer) + "_" + str(ratio) + ".h5", 'a')
            polar_h5file = tables.openFile("hbonds_" + str(isomer) + "_" + str(ratio) + ".h5", 'a')
            for sys in range(10):
                nonpolar_name = isomer + '_' + str(ratio) + '_' + 'inositol_nonpolar' + '_' + str(sys)
                polar_name = isomer + '_' + str(ratio) + '_' + 'inositol_hbonds' + '_' + str(sys)
                
                # get nonpolar matrix
                nonpolar = nonpolar_h5file.root.nonpolar_name
                polar = polar_h5file.root.polar_name
                intersect(nonpolar, polar)
