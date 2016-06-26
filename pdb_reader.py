from math import sqrt




global connectivity
connectivity = {}
#deoxyribose
deoxy = {"C1'": ["C2'", "O4'"],
	"C2'": ["C1'", "C3'", "O2'"],
	"C3'": ["C2'", "C4'", "O3'"],
	"C4'": ["O4'", "C3'", "C5'"],
	"C5'": ["C4'", "O5'"],
	"O3'": ["C3'"],
	"O4'": ["C4'", "C1'"],
	"O5'": ["C5'", "P"],
	"P": ["O5'", 'OP1', 'OP2']}

#cytosine
connectivity['DC'] = {'N1': ['C2', 'C6'], 
	'C6': ['N1', 'C5'],
	'C5': ['C4', 'C6'],
	'C4': ['C5', 'N3', 'N4'],
	'N4': ['C4'],
	'N3': ['C4','C2'],
	'C2': ['O2', 'N1', 'N3'],
	'O2': ['C2']}
#connectivity['DC'] = {**connectivity['DC'], **deoxy}
connectivity['DC'].update(deoxy)
connectivity['DC'].update({'N1': "C1'"})
#adenine
connectivity['DA'] = {'N1': ['C2', 'C6'], 'C6': ['N1', 'N6', 'C5'], 'N6': ['C6'], 'C5': ['C4', 'C6', 'N7'], 'C4': ['N3', 'N9', 'C5'], 'N7': ['C5', 'C8'], 'C8': ['N7', 'N9'], 'N9': ['C4', 'N7']}
connectivity['DA'].update(deoxy)
connectivity['DC'].update({'N9': "C1'"})
#guanine
connectivity['DG'] = {'N1': ['C2', 'C6'], 'C6': ['N1', 'O6', 'C5'], 'O6': ['C6'], 'C5': ['C4', 'C6', 'N7'], 'C4': ['N3', 'N9', 'C5'], 'N7': ['C5', 'C8'], 'C8': ['N7', 'N9'], 'N9': ['C4', 'N7']}
connectivity['DG'].update(deoxy)
connectivity['DC'].update({'N9': "C1'"})
#thymine
connectivity['DT'] = {'N1': ['C2', 'C6'], 'C6': ['N1', 'C5'], 'C5': ['C4', 'C6', 'C7'], 'C4': ['C5', 'O4', 'N3'], 'O4': ['C4'], 'N3': ['C4', 'C2'], 'C2': ['O2', 'N1', 'N3'], 'O2': ['C2']}
connectivity['DT'].update(deoxy)
connectivity['DC'].update({'N1': "C1'"})

class atom:
	def __init__(self, pos, x, y, z, chain='', element='', atom_name='', b=0, charge='', residue='', connections=[]):
		self.position = int(pos)
		self.chain = str(chain)
		self.element = str(element).strip()
		self.atom_name = str(atom_name).strip()
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.b = float(b)
		if charge.strip():
			self.charge = int(charge)
		else:
			self.charge = 0
		self.residue = int(residue)
		self.connections = connections
	def coords(self):
		return(self.x, self.y, self.z)
class residue:
	def __init__(self, index=None, res_name='', chain='', atoms=[]):
		if index:
			self.index = int(index)
		else:
			self.index = None
		self.res_name = str(res_name).strip()
		self.chain = str(chain)
		self.atoms = atoms
		
	def add_atoms(self, atoms):
		self.atoms = atoms
	def get_atoms(self):
		for a in self.atoms:
			yield(a)
		
class structure:
	def __init__(self):
		self.atoms = []
		self.residues = []
		self.connections = {}
	def add_atom(self, atom):
		self.atoms.append(atom)
	def add_residue(self, residue_index, res_name='', chain=''):
		while len(self.residues)< residue_index:
			self.residues.append(residue())
		self.residues[residue_index - 1] = residue(residue_index, res_name, chain)
	def add_atoms_to_residue(self, residue_index, atoms):
		self.residues[residue_index - 1].add_atoms(atoms)
	def connect_residues(res_index1, res_index2):
		pass
	def get_neighbors(distance=5):
		pass
	def residues(self):
		return self.residues
	def get_atoms(self):
		for a in self.atoms:
			yield(a)
	def add_connection(self, atom1, atom2):
		if self.connections.get(atom1.position):
			self.connections[atom1.position].append(atom2.position)
		else:
			self.connections[atom1.position] = [atom2.position]

	def create_connections_by_id(self):
	#create connections by using only the atom IDs
		for residue in self.residues:
			if connectivity.get(residue.res_name):
				for atom1 in residue.atoms:
					for atom2 in residue.atoms:
						if connectivity[residue.res_name].get(atom1.atom_name) and atom2.atom_name in connectivity[residue.res_name][atom1.atom_name]:
							self.add_connection(atom1, atom2)
						elif connectivity[residue.res_name].get(atom2.atom_name) and atom1.atom_name in connectivity[residue.res_name][atom2.atom_name]:
							self.add_connection(atom1, atom2)
	def create_connections_by_distance(self, distance=1.6):
	#create connections by connecting atoms within a certain distance
		a = self.atoms
		for atom1 in range(len(a)):
			for atom2 in range(atom1 + 1, len(a)):
				if sqrt(sum((c1 - c2)**2 for c1, c2 in zip(a[atom1].coords(), a[atom2].coords()))) <= distance:
					self.add_connection(a[atom1], a[atom2])
def read_pdb(pdb_file):
	global connectivity
	pdb = structure()
	residue_no = None
	atoms = []
	for line in open(pdb_file, 'r').readlines():
#ATOM      1  O5'  DC A   1      22.470  37.305  -6.447  1.18 34.97           O
		if line.startswith('ATOM') or line.startswith('HETATOM'):
			if residue_no != int(line[23:26]):
				if residue_no:
					pdb.add_atoms_to_residue(residue_no, atoms)
				atoms = []
				residue_no = int(line[23:26])
				pdb.add_residue(residue_no, line[17:20], line[21:22])
			new_atom = atom(line[6:11], line[30:38], line[38:46], line[46:54], line[21:22], line[76:78], line[12:16], line[60:66], line[78:80], residue_no)
			atoms.append(new_atom)
			pdb.add_atom(new_atom)
	pdb.add_atoms_to_residue(residue_no, atoms)

	return pdb

def unique_dict(d):
	"""
	checks each key and value against each other and removes duplicate associations
	e.g.
	1: [2, 3]
	2: [1, 3, 4, 5]
	becomes
	1: [2, 3]
	2: [4, 5]
	"""
	for k in d.keys():
		for v in d[k]:
			if k in d.get(v):
				d[v].remove(k)
	return d

if __name__ == '__main__':
	pdb = read_pdb('1d89.pdb')
	pdb.create_connections_by_id()
	conns = {}
	for residue in pdb.residues:
		for conn in residue.connections:
		#print(conn, pdb.residues[0].connections[conn])
			conns[conn] = residue.connections[conn]
	#print ("########")
	#print (conns)
	#print ("########")
	#print (unique_dict(conns))
	
	filelist = ['carbon.txt', 'nitrogen.txt', 'oxygen.txt', 'phosphorous.txt', 'all.txt', 'connections.txt']

	f = []
	for filename in filelist:
		f.append(open(filename, 'w'))

	atom_map={}
	atom_map['C'] = 0
	atom_map['N'] = 1
	atom_map['O'] = 2
	atom_map['P'] = 3



	for residue in pdb.atoms:
		index = atom_map.get(residue.element.strip())
		if index != None:
			#f[index].write(str(residue.x) + ';' + str(residue.y) + ';' + str(residue.z) + '\n')
			f[4].write(str(residue.position) + ':' + str(residue.x) + ';' + str(residue.y) + ';' + str(residue.z) + ';' + str(index) + '\n')
		else:
			print(residue.element)



	u = unique_dict(conns)
	for k in u.keys():
		f[5].write(str(k) + ':')
		cells = u[k]
		for x in range(len(cells)):
			cells[x] = str(cells[x])
		f[5].write(';'.join(cells))
		f[5].write('\n')
	for i in range(len(filelist)):
		f[i].close()
