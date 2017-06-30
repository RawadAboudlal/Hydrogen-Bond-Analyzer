'''
Created on May 30, 2017

@author: Rawad
'''

import re
import sys
import time

import loader
from molecule import Bond, HBondType
import numpy as np


# frames should be a dictionary.
def analyzeFrames(frames, maxAngle, maxHBondDistance, maxBondDistance, maxIntermoleculeDistance):
	
	centralAtom1Regex = re.compile("N|O|S|F|C")
	hbondDonors = re.compile("N|O|S|F")
	hbondAcceptors = re.compile("H|N")
	
	# Molecule Id (int): # h-bonds
	totalHBondsCount = {}
	
	# Frame Index (int): list of molecules in longest chain.
	totalHBondChains = {}
	
	# All h-bonds based on which 2 molecules they connect.
	totalHBondsBetweenMolecules = {}
	
	# All h-bonds across all the frames.
	totalHBonds = {}
	
	# frames is a dictionary so we loop through keys this way. key = index of frame; don't have to be in order.
	for frameIndex in frames:
		
		frame = frames[frameIndex]
		
		#molCount = len(frame)
		
		#adjacencyMatrix = []
		
		hbondChains = []
		
		totalHBondChains[frameIndex] = hbondChains
		
		# tuple (sorted, has exactly 2 mol id's): [Bond objects]
		bondsBetweenMolecules = {}
		
		totalHBondsBetweenMolecules[frameIndex] = bondsBetweenMolecules
		
		totalHBonds[frameIndex] = []
		
		for centralMolecule in frame:
			
			if not (centralMolecule.identifier in totalHBondsCount):
				# Each Molecule starts with 0 h-bonds.
				totalHBondsCount[centralMolecule.identifier] = 0
			
			for otherMolecule in frame:
				
				if centralMolecule == otherMolecule:
					continue
				
				withinDistance, otherCentralMoleculeDifference, distanceOtherCentralMolecule = isWithinDistance(centralMolecule.position, otherMolecule.position, maxIntermoleculeDistance)
				
				if not withinDistance:
					continue
				
				del otherCentralMoleculeDifference, distanceOtherCentralMolecule
				
				# Key used for the bonding between these two molecules.
				bondKey = tuple(sorted((centralMolecule.identifier, otherMolecule.identifier)))
				
				# We need to incremenet h-Bond count for both central and other Molecule so this must be initialized.
				if not otherMolecule.identifier in totalHBondsCount:
					totalHBondsCount[otherMolecule.identifier] = 0
				
				for centralAtom1 in centralMolecule.atoms:
					
					# This should be an electronegative atom. EXCEPT for the case of nitrogen being the h-Bond acceptor.
					if not centralAtom1Regex.fullmatch(centralAtom1.element):
						continue
					
					for centralAtom2 in centralMolecule.atoms:
						
						if centralAtom1 == centralAtom2:
							continue
						
						# This is the h-Bond acceptor, needs empty orbital, N can do it too.
						if not hbondAcceptors.fullmatch(centralAtom2.element):
							continue
						
						# Can be bonded across cell boundary.
						withinDistance, central2Central1Diff, distanceCentral2Central1 = isWithinDistance(centralAtom1.position, centralAtom2.position, maxBondDistance)
						
						# Ensures highly electronegative atom and hydrogen are actually bonded.
						if not withinDistance:
							continue
						
						for otherAtom in otherMolecule.atoms:
							
							# This should be an electronegative atom that is being bonded to by the hydrogen.
							if not hbondDonors.fullmatch(otherAtom.element):
								continue
							
							# VMD seems to calculate distance between both electronegative atoms rather than the hydrogen and the other electronegative atom.
							withinDistance, otherCentral1Diff, distanceOtherCentral1 = isWithinDistance(centralAtom1.position, otherAtom.position, maxHBondDistance)
							
							if not withinDistance:
								continue
							
							#print("{} and {} are within distance w/ distance of: {}".format(centralAtom2, otherAtom, distanceOtherCentral2))
							
							cosAngle = np.dot(central2Central1Diff, otherCentral1Diff) / (distanceCentral2Central1 * distanceOtherCentral1)
							angle = np.arccos(cosAngle)
							
							if angle > maxAngle:
								continue
							
							#print("\t{} and {} are also within angle ({}) w/ angle of: {}".format(centralAtom2, otherAtom, maxAngle, angle))
							
							# H-Bond is found here.
							
							# Order of atoms in Bond is important; first one should always be the h-Bond donor (usually hydrogen).
							hbond = Bond(centralAtom2.identifier, otherAtom.identifier, mol1 = centralMolecule, mol2 = otherMolecule)
							
							if not bondKey in bondsBetweenMolecules:
								bondsBetweenMolecules[bondKey] = []
							
							bondsBetweenMolecules[bondKey].append(hbond)
							
							computeChains(hbondChains, hbond)
							
							totalHBondsCount[centralMolecule.identifier] += 1
							totalHBondsCount[otherMolecule.identifier] += 1
							
							totalHBonds[frameIndex].append(hbond)
							
							#print("At frame: {}, Bond: ({}) {}-{} === {} ({})\n\tAs indices: {}-{} === {}".format(frameIndex, centralMolecule.identifier, centralAtom1.element, centralAtom2.element, otherAtom.element, otherMolecule.identifier, centralAtom1.identifier, centralAtom2.identifier, otherAtom.identifier))
							
		
	return (totalHBondsCount, totalHBondChains, totalHBondsBetweenMolecules, totalHBonds)

def outputResult(result, framesCount, outputFileName, maxDistance, maxAngle, hbondTypes, fromFrame, toFrame):
	
	with open(outputFileName, "w") as outputFile:
		
		if toFrame == -1:
			toFrame = fromFrame + framesCount
		
		outputFile.write("From frame {} to frame (exclusive) {} (Total = {})\n".format(fromFrame, toFrame, framesCount))
		
		outputFile.write("Configuration:\n\tMax Hydrogen Bond Distance = {}, Max Hydrogen Bond Angle = {}\n".format(maxDistance, maxAngle))
		
		totalHBondsCount = result[0]
		totalHBondChains = result[1]
		totalHBondsBetweenMolecules = result[2]
		totalHBonds = result[3]
		
		outputFile.write("\n---------- Average HBonds per Molecule ----------\n")
		
		for molId in sorted(totalHBondsCount):
			
			avgHBonds = totalHBondsCount[molId] / framesCount
			
			outputFile.write("Molecule with id {} has an average of {} hydrogen bonds across the simulation.\n".format(molId, avgHBonds))
		
# 		outputFile.write("\n---------- HBonds of Each Molecule ----------\n")
# 		
# 		for frameIndex in totalBondsPerMolecule:
# 			
# 			bondMolDict = totalBondsPerMolecule[frameIndex]
# 			outputFile.write("Frame {}:\n".format(frameIndex))
# 			
# 			for molId in bondMolDict:
# 				outputFile.write("\tMol with id {} has following bonds:\n".format(molId))
# 				outputFile.write("{}\n".format("\n".join("\t\t" + str(Bond) for Bond in bondMolDict[molId])))
		
		outputFile.write("\n---------- HBond Chains ----------\n")
		
		for frameIndex in totalHBondChains:
			
			longestHBondChain = []
			
			for hbondChain in totalHBondChains[frameIndex]:
				
				if len(hbondChain) > len(longestHBondChain):
					longestHBondChain = hbondChain
			
			# TODO: This is temporary.
			if longestHBondChain and len(longestHBondChain) > 3:
				
				idsInChain = []
				
				for hbond in longestHBondChain:
					idsInChain.append(hbond.mol1.identifier)
					idsInChain.append(hbond.mol2.identifier)
				
				uniqueIdsInChain = set(idsInChain)
				
				outputFile.write("Longest hbond chain in frame {}: {}\n".format(frameIndex, "-".join(str(molId) for molId in uniqueIdsInChain)))
				
		
		totalHBondTypes = {}
		
		for hbondType in hbondTypes:
			totalHBondTypes[hbondType.identifier] = 0
		
		for frameIndex in totalHBondsBetweenMolecules:
			
			moleculesToBonds = totalHBondsBetweenMolecules[frameIndex]
			
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
		
		S_HB = {}
		
		initialHBondsList = totalHBonds[fromFrame]
		
		initialHBonds = {i: initialHBondsList[i] for i in range(len(initialHBondsList))}
		
		for initialHBondIndex in initialHBonds:
			S_HB[initialHBondIndex] = 1
		
		for frameIndex in range(fromFrame + 1, toFrame):
			
			currentHBonds = totalHBonds[frameIndex]
			
			# Use new list made from keys to avoid "dictionary changed size during iteration" error.
			for hbondIndex in list(initialHBonds):
				
				hbond = initialHBonds[hbondIndex]
				
				if not isExactHBondPresent(hbond, currentHBonds):
					S_HB[hbondIndex] = 0
					del initialHBonds[hbondIndex]
		
		outputFile.write("\n---------- Continuous Hydrogen Bond Time Correlation Function ----------\n")
		outputFile.write("From frame {} to frame {}\n".format(fromFrame, toFrame))
		
		for hbondIndex in S_HB:
			outputFile.write("HBond {}: {}, S_HB: {}\n".format(hbondIndex, initialHBondsList[hbondIndex], S_HB[hbondIndex]))
		
		outputFile.flush()
	

def isExactHBondPresent(hbondToCheck, hbondChain):
	
	for hbond in hbondChain:
		if hbondToCheck.deep_equal(hbond):
			return True
	
	return False

def computeChains(hbondChains, hbond):

	if hbondChains:
		
		for hbondChain in hbondChains:
			for hbondInChain in hbondChain:
				
				if hbondInChain.is_part_of(hbond.mol1) or hbondInChain.is_part_of(hbond.mol2):
					hbondChain.append(hbond)
					return
			
		
		hbondChains.append([hbond])
		
	else:
		# Make a new chain, with just the one hbond, if there are no chains present.
		hbondChains.append([hbond])
	

# Checks to see if the bonds a Molecule has, bondInMolecule, fit with any of the given types of hbonds, hbondTypes.
def isHBondType(hbondTypes, bondInMolecule):
	
	hbond = HBondType(-1, bondInMolecule)
	
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
				#if x != 0 or y != 0: print("this pos is next to other pos in cell: {}, {}".format(x, y))
				return True, difference, distance
			
			pass
	
	difference = pos2 - pos1
	
	return False, difference, np.linalg.norm(difference)

def main():
	
	try:
		# sys.argv[0] = script name.
		inputFileName = sys.argv[1]
	except:
		print("No input file given.")
		sys.exit(1)
	
	print("Using following input file:", inputFileName)
	
	atomsFileName, atomsPerMolecule, fromFrame, toFrame, maxAngle, maxHBondDistance, maxBondDistance, hbondTypes, outputFileName, maxIntermoleculeDistance = loader.loadInput(inputFileName)
	
	startTime = time.time()
	
	frames = loader.loadAnimatedXyz(atomsFileName, atomsPerMolecule, fromFrame, toFrame)
	
	timeToLoadFrames = time.time() - startTime
	
	minutesToLoadFrames = timeToLoadFrames // 60
	secondsToLoadFrames = timeToLoadFrames % 60
	
	framesCount = len(frames)
	
	print("Loaded {} frames in {minuteMessage}{} seconds.".format(framesCount, secondsToLoadFrames, minuteMessage=(str(minutesToLoadFrames) + " minutes and ") if minutesToLoadFrames > 0 else ""))
	
	startTime = time.time()
	
	result = analyzeFrames(frames, np.radians(maxAngle), maxHBondDistance, maxBondDistance, maxIntermoleculeDistance)
	
	timeToCountHBonds = time.time() - startTime
	
	outputResult(result, framesCount, outputFileName, maxHBondDistance, maxAngle, hbondTypes, fromFrame, toFrame)
	
	minutesToCountHBonds = timeToCountHBonds // 60
	secondsToCountHBonds = timeToCountHBonds % 60
	
	print("Took {minuteMessage}{} seconds to analyze {} frames; that's {} frames/second.".format(secondsToCountHBonds, framesCount, framesCount / timeToCountHBonds, minuteMessage=(str(minutesToCountHBonds) + " minutes and ") if minutesToCountHBonds > 0 else ""))

if __name__ == "__main__":
	main()
