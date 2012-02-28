import os
import numpy
import glob
import logging

bilayer_file = 'popc'
log_file =  'create_systems.log'
LOG_FORMAT = '%(asctime)-15s %(filename)s %(message)s'
logging.basicConfig(filename=log_file, level=logging.DEBUG, format=LOG_FORMAT)
logger = logging.getLogger("mylogger")

def setup(X, Y, Z_bilayer = '15'):
	logger.info("Setting up ... ")
	if not os.path.exists("systems"):
		os.system("mkdir systems")
		logger.info("System directory was created")

	if os.path.exists(log_file):
		# os.system("rm %s" % log_file)
		os.system("rm create_systems.out")
		
	# Create a POPC bilayer with the right XY dimensions if it does not exist
	if not os.path.exists("popc_clean_bounded.gro"):
		command = "editconf_mpi -f popc_clean.gro -o popc_clean_bounded.gro -c -box %(X)s %(Y)s %(Z_bilayer)s" % vars()
		os.system(command)
		logger.info("Bilayer structure not found. Created bilayer file.")
		logger.info(command)
		
		
def cleanup():
	logger.info("Cleaning up")
	os.system("rm systems/\#*")
	os.system("rm systems/temp*.gro")
	os.system("rm \#*")
	logger.info("removed backup files from systems")

def run_log(command):
	logger.info(command)
	os.system(command)

def define_system_geometry(output_dir, num):
	command = "editconf_mpi -f %(output_dir)s/temp.gro -translate 0 0 5.5 -o %(output_dir)s/temp2_%(num)d" % vars()
	run_log(command)

	command = "perl ~/labwork/scripts/analysis/combine_gros.pl %(output_dir)s/temp2_%(num)d.gro popc_clean_bounded.gro %(output_dir)s/test_%(num)d.gro" % vars()
	run_log(command)

	command = "editconf_mpi -f %(output_dir)s/test_%(num)d.gro -o %(output_dir)s/final_%(num)d.gro -translate 0 0 -3" % vars()
	run_log(command)

def add_inositol(isomer, ninos, index, grofile, topfile):
	output = "temp_inositol_%(index)d.gro" % vars()
	os.system("cp system.top systems/%(topfile)s.top" % vars())
	os.system("genbox_mpi -cp systems/%(grofile)s -ci ~/systems/gro/%(isomer)s_em.gro -o systems/%(output)s -p systems/%(topfile)s -nmol %(ninos)d" % vars())
	return output
	
def add_water(index, grofile, topfile):
	os.system("genbox_mpi -cp systems/%(grofile)s -cs ~/systems/gro/tip3.gro -o systems/test_final_%(index)d.gro -p systems/%(topfile)s" % vars())

def main():
	output_dir = 'systems'
	input_dir = 'starting'
	# box dimensions of the POPC bilayer
	popc_dim = ['6.8186', '6.43263', '14']
	X = popc_dim[0]
	Y = popc_dim[1]
	Z = popc_dim[2]
	Z_protein = Z
	setup(X, Y, Z_bilayer = Z)

	num = 0 
	for filename in glob.glob("%s/*.pdb" % input_dir):
		# prepocess the beta oligomer box such that it has the same X and Y dimensions as the POPC bilayer
		newfilename = 'test' + str(num)
		command = "editconf_mpi -f %(filename)s -o %(output_dir)s/temp.gro -box %(X)s %(Y)s %(Z_protein)s -c" % vars()
		run_log(command)
		define_system_geometry(output_dir, num)
		
		grofile = "final_%(num)d.gro" % vars()
		name, ext = os.path.splitext(filename)
		ninos_in_system = int(float(name.split('_')[1]))
		print ninos_in_system
		ninos_add = 32 - ninos_in_system
		temp_file = add_inositol("scyllo", ninos_add, num, grofile, "system_%(num)d" % vars())
		add_water(num, temp_file, "system_%(num)d" % vars())
		num += 1

	cleanup()

if __name__ == '__main__':
	main()
