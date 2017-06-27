'''
Created on May 30, 2017

@author: Rawad
'''

import re

from molecule import atom
from molecule import bond
from molecule import hbond_type
from molecule import molecule

# Will load an xyz file representing an animation. Returns frames as a list; each frame contains a list of molecules.
def loadAnimatedXyz(name, moleculeAtomCount, fromFrame, toFrame):
	
	with open(name) as xyzFile:
		
		atomCountRegex = re.compile(" *(?P<atomCount>[0-9]+) *");
		
		validFloatPattern = "[-+]?[0-9]+[.]*[0-9]*"
		
		atomRegex = re.compile("(?P<element>[A-Za-z]+) *(?P<x>" + validFloatPattern + ") *(?P<y>" + validFloatPattern + ") *(?P<z>" + validFloatPattern + ")")
		
		currentFrame = -1
		
		frames = {}
		
		currentMolecule = molecule(0)
		
		for line in xyzFile:
			
			line = line.strip('\n')
			
			atomCountMatch = atomCountRegex.fullmatch(line)
			
			if atomCountMatch:
				
				currentFrame += 1
				
				if currentFrame < fromFrame:
					#														  + 1 to include the comment line.
					skipLines(xyzFile, int(atomCountMatch.group("atomCount")) + 1)
					# when toFrame = -1, go to end of file (don't skip any lines).
				elif toFrame != -1 and currentFrame >= toFrame:
					return frames
				else:
					
					frames[currentFrame] = []
					
					# Note: this assumes each molecule has the same position in each frame of the file.
					currentMolecule = molecule(0)
				
			
			atomMatch = atomRegex.fullmatch(line)
			
			if atomMatch:
				
				element = atomMatch.group("element")
				
				x = float(atomMatch.group("x"))
				y = float(atomMatch.group("y"))
				z = float(atomMatch.group("z"))
				
				# Makes atom identifier 1-indexed. Also assumes each atom has the same position in each molecule.
				a = atom(element, len(currentMolecule.atoms) + 1, x, y, z);
				
				currentMolecule.atoms.append(a)
				
				if len(currentMolecule.atoms) >= moleculeAtomCount:
					
					# Calculates approximate center-of-mass of molecule once it is fully loaded.
					currentMolecule.calculateCenter()
					
					frames[currentFrame].append(currentMolecule)
					currentMolecule = molecule(len(frames[currentFrame]))
				
		return frames
	
	print("Couldn't open file", name)
	
	return None

def skipLines(file, linesToSkip):
	
	for i in range(linesToSkip):
		file.readline()
		
		i = i# Annoying unused variable warning.
	
def loadInput(name):
	
	with open(name) as inputFile:
		
		atomsFileRegex = re.compile("atoms file", re.IGNORECASE)
		atomsPerMoleculeRegex = re.compile("atoms per molecule", re.IGNORECASE)
		fromFrameRegex = re.compile("from frame", re.IGNORECASE)
		toFrameRegex = re.compile("to frame", re.IGNORECASE)
		maxAngleRegex = re.compile("max angle", re.IGNORECASE)
		maxHBondDistanceRegex = re.compile("max h-bond distance", re.IGNORECASE)
		maxBondDistanceRegex = re.compile("max bond distance", re.IGNORECASE)
		hbondTypesRegex = re.compile("hbond types", re.IGNORECASE)
		outputFileRegex = re.compile("output file", re.IGNORECASE)
		
		# .+? -> ? matches up to the FIRST closing square bracker ].
		listRegex = re.compile("(?P<list>\[.+?\])")
		bondRegex = re.compile("(?P<atom1>[0-9]+)-(?P<atom2>[0-9]+)")
		
		atomsFile = ""
		atomsPerMolecule = 1
		fromFrame = 0
		toFrame = -1
		maxAngle = 0
		maxHBondDistance = 1
		maxBondDistance = 1
		hbondTypes = []
		outputFileName = ""
		
		for line in inputFile:
			
			if line.startswith("#"):
				continue
			
			try:
				
				# Ignored comments at end of line if present.
				tokens = line.partition("#")[0].split("=")
				
				key = tokens[0].strip()
				value = tokens[1].strip()
				
			except:
				print("Error occured parsing line:", line, end="")
				continue
			
			if atomsFileRegex.fullmatch(key):
				atomsFile = str(value)
			elif atomsPerMoleculeRegex.fullmatch(key):
				atomsPerMolecule = int(value)
			elif fromFrameRegex.fullmatch(key):
				fromFrame = int(value)
			elif toFrameRegex.fullmatch(key):
				toFrame = int(value)
			elif maxAngleRegex.fullmatch(key):
				maxAngle = float(value)
			elif maxHBondDistanceRegex.fullmatch(key):
				maxHBondDistance = float(value)
			elif maxBondDistanceRegex.fullmatch(key):
				maxBondDistance = float(value)
			elif hbondTypesRegex.fullmatch(key):
				
				for lstMatch in listRegex.finditer(value):
					
					bonds = []
					
					for bondMatch in bondRegex.finditer(lstMatch.group("list")):
						bonds.append(bond(int(bondMatch.group("atom1")), int(bondMatch.group("atom2"))))
					
					# Makes hbond_type identifiers 1-indexed.
					hbondType = hbond_type(len(hbondTypes) + 1, bonds)
					hbondTypes.append(hbondType)
				
			elif outputFileRegex.fullmatch(key):
				outputFileName = value
					
		return atomsFile, atomsPerMolecule, fromFrame, toFrame, maxAngle, maxHBondDistance, maxBondDistance, hbondTypes, outputFileName
