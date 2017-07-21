'''
Created on May 30, 2017

@author: Rawad
'''

import re
import sys
import time
import os

import loader
from loop import Graph
from molecule import Bond, HBondType
import numpy as np


# frames should be a dictionary.
def analyzeFrames(frames, maxAngle, maxHBondDistance, maxBondDistance, maxIntermoleculeDistance):
	
	centralAtom1Regex = re.compile("N|O|S|F")# "N|O|S|F|C". Keep this way for now.
	hbondDonors = re.compile("N|O|S|F")
	hbondAcceptors = re.compile("H|N")
	
	# Molecule Id (int): # h-bonds
	totalHBondsCount = {}
	
	# Frame index (int): Graph
	hbondGraphs = {}
	
	# All h-bonds based on which 2 molecules they connect.
	totalHBondsBetweenMolecules = {}
	
	# All h-bonds across all the frames.
	totalHBonds = {}
	
	# frames is a dictionary so we loop through keys this way. key = index of frame; don't have to be in order.
	for frameIndex in frames:
		
		frame = frames[frameIndex]
		
		graph = Graph()
		
		hbondGraphs[frameIndex] = graph
		
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
							
							graph.addGroup(hbond.mol1.identifier, hbond.mol2.identifier)
							
							totalHBondsCount[centralMolecule.identifier] += 1
							totalHBondsCount[otherMolecule.identifier] += 1
							
							totalHBonds[frameIndex].append(hbond)
							
							#print("At frame: {}, Bond: ({}) {}-{} === {} ({})\n\tAs indices: {}-{} === {}".format(frameIndex, centralMolecule.identifier, centralAtom1.element, centralAtom2.element, otherAtom.element, otherMolecule.identifier, centralAtom1.identifier, centralAtom2.identifier, otherAtom.identifier))
							
		
	return (totalHBondsCount, hbondGraphs, totalHBondsBetweenMolecules, totalHBonds)

def outputResult(result, framesCount, maxDistance, maxAngle, hbondTypes, fromFrame, toFrame):
	
	totalHBondsCount = result[0]
	hbondGraphs = result[1]
	totalHBondsBetweenMolecules = result[2]
	totalHBonds = result[3]
	
	with open("hbonds.txt", "w") as outputFile:
		
		if toFrame == -1:
			toFrame = fromFrame + framesCount
		
		outputFile.write("From frame {} to frame {} (Total = {})\n".format(fromFrame, toFrame, framesCount))
		
		outputFile.write("Configuration:\n\tMax Hydrogen Bond Distance = {}, Max Hydrogen Bond Angle = {}\n".format(maxDistance, maxAngle))
		
		outputFile.write("\n---------- Average HBonds per Molecule ----------\n")
		
		for molId in sorted(totalHBondsCount):
			
			avgHBonds = totalHBondsCount[molId] / framesCount
			
			outputFile.write("Molecule with id {} has an average of {:.2f} hydrogen bonds across the simulation.\n".format(molId, avgHBonds))
		
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
		
		connectedHBondComponents = {}
		
		for frameIndex in hbondGraphs:
			connectedHBondComponents[frameIndex] = hbondGraphs[frameIndex].getConnectedComponents()
		
		outputFile.write("\n---------- HBond Chains ----------\n")
		
		hbondChainsCount = {}
		
		for frameIndex in connectedHBondComponents:
			
			for hbondChain in connectedHBondComponents[frameIndex]:
				
				moleculesInChain = tuple(mol for mol in hbondChain.graph)
				
				if moleculesInChain in hbondChainsCount:
					hbondChainsCount[moleculesInChain] += 1
				else:
					hbondChainsCount[moleculesInChain] = 1
		
		numberOfUniqueChains = len(hbondChainsCount)
		
		if numberOfUniqueChains != 0:
			numberOfTotalChains = 0
			
			for hbondChain in hbondChainsCount:
				numberOfTotalChains += hbondChainsCount[hbondChain] * len(hbondChain)
			
			outputFile.write("Average chains length: {:.2f}\n".format(numberOfTotalChains / numberOfUniqueChains))	
			
			outputFile.write("\n---------- HBond Loops ----------\n")
			
			hbondLoopsCount = {}
			
			for frameIndex in connectedHBondComponents:
				
				hbondChains = connectedHBondComponents[frameIndex]
				
				for hbondChain in hbondChains:
					if hbondChain.isCyclic():
						
						hbondLoop = tuple(mol for mol in hbondChain.graph)
						
						if hbondLoop in hbondLoopsCount:
							hbondLoopsCount[hbondLoop] += 1
						else:
							hbondLoopsCount[hbondLoop] = 1
						
		
		numberOfUniqueLoops = len(hbondLoopsCount)
		
		if numberOfUniqueLoops != 0:
			
			numberOfTotalLoops = 0
			
			for hbondLoop in hbondLoopsCount:
				numberOfTotalLoops += hbondLoopsCount[hbondLoop] * len(hbondLoop)
			
			outputFile.write("The average size of hbond loops is: {:.2f}\n".format(numberOfTotalLoops / numberOfUniqueLoops))
		
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
		
		outputFile.flush()
		
	with open("lifetime.dat", "w") as outputFile:
		
		# Split frames into 5 time origins.
		segments = 5
		
		outputFile.write("#---------- Continuous Hydrogen Bond Time Correlation Function ----------\n")
		outputFile.write("#From frame {} to frame {}. Number of origins: {}.\n".format(fromFrame, toFrame, segments))
		outputFile.write("#Frame Index | S_HB\n")
		
		S_HB_sum = {}
		
		segmentSize = framesCount // segments
		
		for segmentIndex in range(segments):
			
			newStartFrame = segmentSize * segmentIndex + fromFrame
			
			initialHBondsList = totalHBonds[newStartFrame]
			
			h_of_0 = len(initialHBondsList)
			
			H_of_t = {}
			
			initialHBonds = {i: initialHBondsList[i] for i in range(len(initialHBondsList))}
			
			for initialHBondIndex in initialHBonds:
				H_of_t[initialHBondIndex] = 1
			
			# Don't add +1 to end. e.g. 500 frames (0-499), will go from 0-99, 100-199, 200-299, 300-399, and 400-499.
			for t in range(newStartFrame, newStartFrame + segmentSize):
				
				currentHBonds = totalHBonds[t]
				
				# Use new list made from keys to avoid "dictionary changed size during iteration" error.
				for hbondIndex in list(initialHBonds):
					
					hbond = initialHBonds[hbondIndex]
					
					if not isExactHBondPresent(hbond, currentHBonds):
						H_of_t[hbondIndex] = 0
						del initialHBonds[hbondIndex]
				
				# Make it so final index is relative to 0; i.e. the first t will always be 0.
				S_HB_index = t - newStartFrame
				
				if not S_HB_index in S_HB_sum:
					S_HB_sum[S_HB_index] = 0
				
				S_HB_sum[S_HB_index] += sum(H_of_t[hbondIndex] for hbondIndex in H_of_t) / h_of_0
		
		for t in S_HB_sum:
			# / segments to average out S_HB function over the given
			outputFile.write("{} {}\n".format(t, S_HB_sum[t] / segments))
		
		outputFile.flush()
	
	with open("chain_size.dat", "w") as outputFile:
		
		outputFile.write("#---------- Chain Size Frequency ----------\n")
		outputFile.write("#From frame {} to frame {}.\n".format(fromFrame, toFrame))
		outputFile.write("#Chain Size | Number of Occurences\n")
		
		hbondChainSizes = {}
		
		for hbondChain in hbondChainsCount:
			hbondChainSize = len(hbondChain)
			
			if hbondChainSize in hbondChainSizes:
				hbondChainSizes[hbondChainSize] += hbondChainsCount[hbondChain]
			else:
				hbondChainSizes[hbondChainSize] = hbondChainsCount[hbondChain]
			
		
		for hbondChainSize in hbondChainSizes:
			outputFile.write("{} {}\n".format(hbondChainSize, hbondChainSizes[hbondChainSize]))
		
		outputFile.flush()
	

def isExactHBondPresent(hbondToCheck, hbondChain):
	
	for hbond in hbondChain:
		if hbondToCheck.deep_equal(hbond):
			return True
	
	return False

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
	
	atomsFileName, atomsPerMolecule, fromFrame, toFrame, maxAngle, maxHBondDistance, maxBondDistance, hbondTypes, maxIntermoleculeDistance = loader.loadInput(inputFileName)
	
	startTime = time.time()
	
	frames = loader.loadAnimatedXyz(atomsFileName, atomsPerMolecule, fromFrame, toFrame)
	
	timeToLoadFrames = time.time() - startTime
	
	minutesToLoadFrames = timeToLoadFrames // 60
	secondsToLoadFrames = timeToLoadFrames % 60
	
	framesCount = len(frames)
	
	print("Loaded {} frames in {minuteMessage}{:.2f} seconds.".format(framesCount, secondsToLoadFrames, minuteMessage=("{} minutes and ".format(minutesToLoadFrames)) if minutesToLoadFrames > 0 else ""))
	
	startTime = time.time()
	
	result = analyzeFrames(frames, np.radians(maxAngle), maxHBondDistance, maxBondDistance, maxIntermoleculeDistance)
	
	# Extract directory from input file and use it to save output files. Doing it down here won't affect loading the input file.
	os.chdir(inputFileName.rsplit('/', 1)[0])
	
	outputResult(result, framesCount, maxHBondDistance, maxAngle, hbondTypes, fromFrame, toFrame)
	
	timeToCountHBonds = time.time() - startTime
	
	minutesToCountHBonds = timeToCountHBonds // 60
	secondsToCountHBonds = timeToCountHBonds % 60
	
	print("Took {minuteMessage}{:.2f} seconds to analyze {} frames; that's {:.2f} frames/second.".format(secondsToCountHBonds, framesCount, framesCount / timeToCountHBonds, minuteMessage=("{} minutes and ".format(minutesToCountHBonds)) if minutesToCountHBonds > 0 else ""))

if __name__ == "__main__":
	main()
