'''
Created on May 30, 2017

@author: Rawad
'''

import numpy as np

class Molecule:
	
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
	
class Atom:
	
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
	

class Bond:
	
	# Can represent h-Bond between 2 molecules or a template Bond for h-bonds between generic molecules or a covalent Bond within the same Molecule.
	def __init__(self, atom1, atom2, mol1 = None, mol2 = None):
		self.atom1 = atom1
		self.atom2 = atom2
		self.mol1 = mol1
		self.mol2 = mol2
	
	def is_part_of(self, molecule):
		'''
		Returns True if molecule is the same as one of the molecules (mol1 or mol1) that make up this Bond.
		'''
		return self.mol1.identifier == molecule.identifier or self.mol2.identifier == molecule.identifier
	
	def deep_equal(self, other):
		'''
		Must satisfy self.__eq__(other) and must also have the same molecules bonded as other.
		'''
		return self.__eq__(other) and sorted((self.mol1.identifier, self.mol2.identifier)) == sorted((other.mol1.identifier, other.mol2.identifier))
	
	def __eq__(self, other):
		'''
		Important Note: Assumes that the atoms are sorted such that atom1 is the H of the h-Bond and atom2 is the
		electronegative Atom the H is bonded with.
		'''
		return self.atom1 == other.atom1 and self.atom2 == other.atom2
	
	def __ne__(self, other):
		return not self.__eq__(other)
	
	def __repr__(self):
		return "Bond(Id 1: {}, Id 2: {})".format(self.atom1, self.atom2)
	

class HBondType:
	
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
		
		for Bond in other.bonds:
			if Bond in self.bonds:
				matches += 1
			else:
				return False
		
		return matches >= numBonds
		
	
	def __repr__(self):
		return "HBondType(Bonds: {})".format(self.bonds)
	
