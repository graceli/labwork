import re
import sys


class System:
	def __init__(self):
		self.comment = ""
		self.num_atoms = ""
		self.box_size = ""
		self.protein = []
		self.solvent = []
		self.solute = []
		self.lipid = []
	
	def to_string(self):
		""" returns in the gro file format in the order of protein, lipid, ins, solvent"""	
		protein = '\n'.join(self.protein)
		lipid = '\n'.join(self.lipid)
		solute = '\n'.join(self.solute)
		solvent = '\n'.join(self.solvent)
		box_size = self.box_size
		contents = self.comment + '\n' + str(self.num_atoms) + '\n' + protein + '\n' + solute + '\n' + lipid + '\n' + solvent + '\n' + box_size
		return contents	
	
	def to_string_no_solute(self):
		protein = '\n'.join(self.protein)
		lipid = '\n'.join(self.lipid)
		solvent = '\n'.join(self.solvent)
		box_size = self.box_size
		num_atoms = len(self.protein) + len(self.lipid) + len(self.solvent)
		contents = self.comment + '\n' + str(num_atoms) + '\n' + protein + '\n' + lipid + '\n' + solvent + '\n' + box_size
		return contents	
	

		
# process the gro file
system = System()
num_line = 0
for line in sys.stdin:
	line = line.rstrip()
	if num_line == 0:
		system.comment = line
		num_line += 1
		continue
	
	if num_line == 1:
		system.num_atoms = line	
		num_line += 1
		continue
	
	if num_line == int(system.num_atoms) + 2:
		system.box_size = line
		break
	
	fields = line.split()	
	resname = fields[0]
	if re.match(".*SOL.*", resname): 
		system.solvent.append(line)
	elif re.match(".*POP.*", resname): 
		system.lipid.append(line)
	elif re.match(".*INS.*", resname):
		system.solute.append(line)
	else:
		system.protein.append(line)

	num_line += 1

#print system.to_string()	
print system.to_string_no_solute()
