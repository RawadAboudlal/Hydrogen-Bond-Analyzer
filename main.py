'''
Created on May 30, 2017

@author: Rawad
'''

import re
import sys
import time

import loader
from molecule import bond, hbond_type
import numpy as np


# frames should be a dictionary.
def analyzeFrames(frames, maxAngle, maxHBondDistance, maxBondDistance):
	
	centralAtom1Regex = re.compile("N|O|S|F|C")
	hbondDonors = re.compile("N|O|S|F")
	hbondAcceptors = re.compile("H|N")
	
	# Molecule Id (int): # h-bonds
	totalHBondsCount = {}
	
	# Frame Index (int): list of molecules in longest chain.
	hbondChains = {}
	
	# All bonds, for each molecule, across all the frames.
	totalHBonds = {}
	
	# frames is a dictionary so we loop through keys this way. key = index of frame; don't have to be in order.
	for frameIndex in frames:
		
		frame = frames[frameIndex]
		
		# Holds list of all the chains
		chainsList = []
		
		molCount = len(frame)
		
		adjacencyMatrix = []
		
		# tuple (sorted, has exactly 2 mol id's): [bond objects]
		bondsBetweenMolecules = {}
		
		totalHBonds[frameIndex] = bondsBetweenMolecules
		
		for centralMolecule in frame:
			
			if not (centralMolecule.identifier in totalHBondsCount):
				# Each molecule starts with 0 h-bonds.
				totalHBondsCount[centralMolecule.identifier] = 0
			
			for centralAtom1Index in range(len(centralMolecule.atoms)):
				
				centralAtom1 = centralMolecule.atoms[centralAtom1Index]
				
				# This should be an electronegative atom. EXCEPT for the case of nitrogen being the h-bond acceptor.
				if not centralAtom1Regex.fullmatch(centralAtom1.element):
					continue
				
				# +1 so we don't try to match this atom with itself.
				for centralAtom2Index in range(len(centralMolecule.atoms)):
					
					centralAtom2 = centralMolecule.atoms[centralAtom2Index]
					
					if centralAtom1 == centralAtom2:
						continue
					
					# This is the h-bond acceptor, needs empty orbital, N can do it too.
					if not hbondAcceptors.fullmatch(centralAtom2.element):
						continue
					
					# Can be bonded across cell boundary.
					withinDistance, central2Central1Diff, distanceCentral2Central1 = isWithinDistance(centralAtom1.position, centralAtom2.position, maxBondDistance)
					
					# Ensures highly electronegative atom and hydrogen are actually bonded.
					if not withinDistance:
						continue
					
					# +1 so we don't match molecule with itself.
					for otherMolecule in frame:
						
						if centralMolecule == otherMolecule:
							continue
						
						# Key used for the bonding between these two molecules.
						bondKey = tuple(sorted((centralMolecule.identifier, otherMolecule.identifier)))
						
						# We need to incremenet h-bond count for both central and other molecule so this must be initialized.
						if not otherMolecule.identifier in totalHBondsCount:
							totalHBondsCount[otherMolecule.identifier] = 0
						
						for otherAtom in otherMolecule.atoms:
							
							# This should be an electronegative atom that is being bonded to by the hydrogen.
							if not hbondDonors.fullmatch(otherAtom.element):
								continue
							
							withinDistance, otherCentral2Diff, distanceOtherCentral2 = isWithinDistance(centralAtom2.position, otherAtom.position, maxHBondDistance)
							
							if not withinDistance:
								continue
							
							cosAngle = np.dot(central2Central1Diff, otherCentral2Diff) / (distanceCentral2Central1 * distanceOtherCentral2)
							angle = np.arccos(cosAngle)
							
							if angle > maxAngle:
								continue
							
							# H-bond is found here.
							
							# Order of atoms in bond is important; first one should always be the h-bond acceptor.
							b = bond(centralAtom2.identifier, otherAtom.identifier, mol1 = centralMolecule, mol2 = otherMolecule)
							
							if not bondKey in bondsBetweenMolecules:
								bondsBetweenMolecules[bondKey] = []
							
							bondsBetweenMolecules[bondKey].append(b)
							
							for hbondChain in hbondChains[frameIndex]:
								
								if centralMolecule.identifier in hbondChain:
									hbondChain.append(otherMolecule.identifier)
									break
								elif otherMolecule.identifier in hbondChain:
									hbondChain.append(centralMolecule.identifier)
									break
								
							else:
								hbondChains[frameIndex].append([centralMolecule.identifier, otherMolecule.identifier])
							
							totalHBondsCount[centralMolecule.identifier] += 1
							totalHBondsCount[otherMolecule.identifier] += 1
							
		
	return (totalHBondsCount, hbondChains, totalHBonds)

def outputResult(result, framesCount, outputFileName, maxDistance, maxAngle, hbondTypes):
	
	with open(outputFileName, "w") as outputFile:
		
		outputFile.write("Total frames: {}\n".format(framesCount))
		
		outputFile.write("Configuration:\n\tMax Hydrogen Bond Distance = {}, Max Hydrogen Bond Angle = {}\n".format(maxDistance, maxAngle))
		
		totalHBondsCount = result[0]
		hbondChains = result[1]
		totalHBonds = result[2]
		
		outputFile.write("\n---------- Average HBonds per Molecule ----------\n")
		
		for molId in totalHBondsCount:
			
			avgHBonds = totalHBondsCount[molId] / framesCount
			
			outputFile.write("Molecule with id {} has an average of {} hydrogen bonds across the simulation.\n".format(molId, avgHBonds))
		
# 		outputFile.write("\n---------- HBonds Chains ----------\n")
# 		
# 		for frameIndex in hbondChains:
# 			outputFile.write("Frame {}:\n".format(frameIndex))
# 			for hbondChain in hbondChains[frameIndex]:
# 				outputFile.write("\tMolecules in this Chain: {}\n".format(", ".join(str(molId) for molId in set(hbondChain))))
		
# 		outputFile.write("\n---------- HBonds of Each Molecule ----------\n")
# 		
# 		for frameIndex in totalBondsPerMolecule:
# 			
# 			bondMolDict = totalBondsPerMolecule[frameIndex]
# 			outputFile.write("Frame {}:\n".format(frameIndex))
# 			
# 			for molId in bondMolDict:
# 				outputFile.write("\tMol with id {} has following bonds:\n".format(molId))
# 				outputFile.write("{}\n".format("\n".join("\t\t" + str(bond) for bond in bondMolDict[molId])))
		
		totalHBondTypes = {}
		
		for hbondType in hbondTypes:
			totalHBondTypes[hbondType.identifier] = 0
		
		for frameIndex in totalHBonds:
			
			moleculesToBonds = totalHBonds[frameIndex]
			
			for moleculePair in moleculesToBonds:
				
				bonds = moleculesToBonds[moleculePair]
				
				#print("{} have {} bonds with each other.".format(moleculePair, len(bonds)))
				
				hbondTypeMatch = isHBondType(hbondTypes, bonds)
				
				if not hbondTypeMatch == None:
					#print("Found hbond with type {} and mol id {}".format(hbondTypeMatch.identifier, molId))
					totalHBondTypes[hbondTypeMatch.identifier] += 1
				
		
		outputFile.write("\n---------- Total Number of HBonds ----------\n")
		
		for hbondType in totalHBondTypes:
			outputFile.write("HBond Type {} was found {} time(s).\n".format(hbondType, totalHBondTypes[hbondType]))
		
		
		outputFile.flush()
	

# Checks to see if the bonds a molecule has, bondInMolecule, fit with any of the given types of hbonds, hbondTypes.
def isHBondType(hbondTypes, bondInMolecule):
	
	hbond = hbond_type(-1, bondInMolecule)
	
	for hbondType in hbondTypes:
		if hbondType.matches(hbond):
			return hbondType
	
	return None

# This will return the following values: boolean, pos2 - pos1, distance. First is whether two positions are within the given
# maxDistance. Second is the difference between pos2 - pos1 OR pos2 in the image cell that is within maxDistance of pos1.
# Third is the distance of pos2 - pos1.
def isWithinDistance(pos1, pos2, maxDistance, La=np.array([34.60467452, 0, 0]), Lb=np.array([17.30233726, 29.968527222, 0])):
	
	# This can instead be loaded from input file.
	#[34.60467452, 0, 0],
	#[17.30233726, 29.968527222, 0],
	#[0, 0, 100] # This is unused because we are only interested in neighbouring cells in the xy plane.
	
	cellsToSearch = (-1, 0, 1)
	
	for x in cellsToSearch:
		for y in cellsToSearch:
			
			newPos2 = pos2 + (x * La) + (y * Lb)# + (z * Lc)
			
			difference = newPos2 - pos1
			distance = np.linalg.norm(difference)
			
			if distance <= maxDistance:
				return True, difference, distance
			
			pass
	
	difference = pos2 - pos1
	
	return False, difference, np.linalg.norm(difference)

def main():
	
	try:
		# sys.argv[0] = script name.
		inputFileName = sys.argv[1]
	except:
		print("Not input file given.")
		sys.exit(1)
	
	print("Using following input file:", inputFileName)
	
	atomsFileName, atomsPerMolecule, fromFrame, toFrame, maxAngle, maxHBondDistance, maxBondDistance, hbondTypes, outputFileName = loader.loadInput(inputFileName)
	
	startTime = time.time()
	
	frames = loader.loadAnimatedXyz(atomsFileName, atomsPerMolecule, fromFrame=fromFrame, toFrame=toFrame)
	
	timeToLoadFrames = time.time() - startTime
	
	minutesToLoadFrames = timeToLoadFrames // 60
	secondsToLoadFrames = timeToLoadFrames % 60
	
	framesCount = len(frames)
	
	print("Loaded {} frames in {minuteMessage}{} seconds.".format(framesCount, secondsToLoadFrames, minuteMessage=(str(minutesToLoadFrames) + " minutes and ") if minutesToLoadFrames > 0 else ""))
	
	startTime = time.time()
	
	result = analyzeFrames(frames, np.radians(maxAngle), maxHBondDistance, maxBondDistance)
	
	timeToCountHBonds = time.time() - startTime
	
	outputResult(result, framesCount, outputFileName, maxHBondDistance, maxAngle, hbondTypes)
	
	minutesToCountHBonds = timeToCountHBonds // 60
	secondsToCountHBonds = timeToCountHBonds % 60
	
	print("Took {minuteMessage}{} seconds to analyze {} frames; that's {} frames/second.".format(secondsToCountHBonds, framesCount, framesCount / timeToCountHBonds, minuteMessage=(str(minutesToCountHBonds) + " minutes and ") if minutesToCountHBonds > 0 else ""))

if __name__ == "__main__":
	main()
