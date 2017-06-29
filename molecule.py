'''
Created on May 30, 2017

@author: Rawad
'''

import numpy as np

class molecule:
	
	def __init__(self, identifier):
		self.atoms = []
		self.identifier = identifier
		
	def __repr__(self):
		return "Molecule({})".format(self.identifier)
	
	def calculateCenter(self):
		
		atomsCount = len(self.atoms)
		
		self.x = sum(a.x for a in self.atoms) / atomsCount
		self.y = sum(a.y for a in self.atoms) / atomsCount
		self.z = sum(a.z for a in self.atoms) / atomsCount
		
		self.position = np.array([self.x, self.y, self.z])
	
class atom:
	
	def __init__(self, element, identifier, x, y, z):
		self.element = element
		self.identifier = identifier
		self.x = x
		self.y = y
		self.z = z
		self.position = np.array([x, y, z])
	
	def __eq__(self, other):
		return self.identifier == other.identifier
	
	def __ne__(self, other):
		return not self.__eq__(other)
	
	def __repr__(self):
		return "Atom({}, {}, {}, {}, {})".format(self.element, self.identifier, self.x, self.y, self.z)
	

class bond:
	
	# Can represent h-bond between 2 molecules or a template bond for h-bonds between generic molecules or a covalent bond within the same molecule.
	def __init__(self, atom1, atom2, mol1 = None, mol2 = None):
		self.atom1 = atom1
		self.atom2 = atom2
		self.mol1 = mol1
		self.mol2 = mol2
	
	def is_part_of(self, molecule):
		return self.mol1.identifier == molecule.identifier or self.mol2.identifier == molecule.identifier
	
	def __eq__(self, other):
		'''
		Important Note: Assumes that the atoms are sorted such that atom1 is the H of the h-bond and atom2 is the
		electronegative atom the H is bonded with.
		'''
		return self.atom1 == other.atom1 and self.atom2 == other.atom2
	
	def __ne__(self, other):
		return not self.__eq__(other)
	
	def __repr__(self):
		return "Bond(Id 1: {}, Id 2: {})".format(self.atom1, self.atom2)
	

class hbond_type:
	
	# bonds = list of bonds that characterize this hbond between 2 molecules.
	def __init__(self, identifier, bonds):
		self.identifier = identifier
		self.bonds = bonds
	
	def matches(self, other):
		'''
		Checks to see if some combinations of bonds in other matches ALL bonds in self. Always call this from the template
		type so that other is the one obtained from analyzing the simulation.
		'''
		
		numBonds = len(self.bonds)
		
		matches = 0
		
		for bond in other.bonds:
			if bond in self.bonds:
				matches += 1
			else:
				return False
		
		return matches >= numBonds
		
	
	def __repr__(self):
		return "HBondType(Bonds: {})".format(self.bonds)
	
