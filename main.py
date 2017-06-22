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
def analyzeFrames(frames, maxAngle, maxHBondDistance, maxBondDistance, hbondTypes):
	
	electronegativeAtoms = re.compile("N|O|S|F")
	
	# Molecule Id (int): # h-bonds
	totalHBonds = {}
	
	totalHBondTypes = {}
	
	for hbondType in hbondTypes:
		totalHBondTypes[hbondType.identifier] = 0
	
	# Frame Index (int): list of molecules in longest chain.
	hbondChains = {}
	
	# All bonds, for each molecule, across all the frames.
	totalHBondsPerMolecule = {}
	
	# frames is a dictionary so we loop through keys this way. key = index of frame; don't have to be in order.
	for frameIndex in frames:
		
		frame = frames[frameIndex]
		
		hbondChains[frameIndex] = []
		
		molCount = len(frame)
		
		adjacencyMatrix = []
		
		bondsPerMolecule = {}
		
		totalHBondsPerMolecule[frameIndex] = bondsPerMolecule
		
		for centralMolecule in frame:
			
			if not (centralMolecule.identifier in totalHBonds):
				# Each molecule starts with 0 h-bonds.
				totalHBonds[centralMolecule.identifier] = 0
			
			bondsPerMolecule[centralMolecule.identifier] = []
			
			for centralAtom1Index in range(len(centralMolecule.atoms)):
				
				centralAtom1 = centralMolecule.atoms[centralAtom1Index]
				
				# This should be an electronegative atom.
				if not electronegativeAtoms.fullmatch(centralAtom1.element):
					continue
				
				# +1 so we don't try to match this atom with itself.
				for centralAtom2Index in range(len(centralMolecule.atoms)):
					
					centralAtom2 = centralMolecule.atoms[centralAtom2Index]
					
					if centralAtom1 == centralAtom2:
						continue
					
					# This should always be a hydrogen
					if not centralAtom2.element == "H":
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
						
						if not otherMolecule.identifier in bondsPerMolecule:
							bondsPerMolecule[otherMolecule.identifier] = []
						
						# We need to incremenet h-bond count for both central and other molecule so this must be initialized.
						if not otherMolecule.identifier in totalHBonds:
							totalHBonds[otherMolecule.identifier] = 0
						
						for otherAtom in otherMolecule.atoms:
							
							# This should be an electronegative atom that is being bonded to by the hydrogen.
							if not electronegativeAtoms.fullmatch(otherAtom.element):
								continue
							
							withinDistance, otherCentral2Diff, distanceOtherCentral2 = isWithinDistance(centralAtom2.position, otherAtom.position, maxHBondDistance)
							
							if not withinDistance:
								continue
							
							cosAngle = np.dot(central2Central1Diff, otherCentral2Diff) / (distanceCentral2Central1 * distanceOtherCentral2)
							angle = np.arccos(cosAngle)
							
							if angle > maxAngle:
								continue
							
							# H-bond is found here.
							
							b = bond(centralAtom2.identifier, otherAtom.identifier, mol1 = centralMolecule, mol2 = otherMolecule)
							
							bondsPerMolecule[centralMolecule.identifier].append(b)
							bondsPerMolecule[otherMolecule.identifier].append(b)
							
							for hbondChain in hbondChains[frameIndex]:
								
								if centralMolecule.identifier in hbondChain:
									hbondChain.append(otherMolecule.identifier)
									break
								elif otherMolecule.identifier in hbondChain:
									hbondChain.append(centralMolecule.identifier)
									break
								
							else:
								hbondChains[frameIndex].append([centralMolecule.identifier, otherMolecule.identifier])
							
							totalHBonds[centralMolecule.identifier] += 1
							totalHBonds[otherMolecule.identifier] += 1
							
		
		# At end of each frame, go through all the found h-bonds and test them against the given h-bond types
		for molId in bondsPerMolecule:
			
			hbond = hbond_type(-1, bondsPerMolecule[molId])
			
			for hbondType in hbondTypes:
				if hbondType.matches(hbond):
					totalHBondTypes[hbondType.identifier] += 1
		
	return (totalHBonds, hbondChains, totalHBondsPerMolecule, totalHBondTypes)

# This will return the following values: boolean, pos2 - pos1, distance. First is whether two positions are within the given
# maxDistance. Second is the difference between pos2 - pos1 OR pos2 in the image cell that is within maxDistance of pos1.
# Third is the distance of pos2 - pos1.
def isWithinDistance(pos1, pos2, maxDistance, La=np.array([34.60467452, 0, 0]), Lb=np.array([17.30233726, 29.968527222, 0])):
	
	#[34.60467452, 0, 0],
	#[17.30233726, 29.968527222, 0],
	#[0, 0, 100]
	
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

def outputResult(result, framesCount, outputFileName, maxDistance, maxAngle):
	
	with open(outputFileName, "w") as outputFile:
		
		outputFile.write("Total frames: {}\n".format(framesCount))
		
		outputFile.write("Configuration:\n\tMax Hydrogen Bond Distance = {}, Max Hydrogen Bond Angle = {}\n".format(maxDistance, maxAngle))
		
		totalHBonds = result[0]
		hbondChains = result[1]
		totalBondsPerMolecule = result[2]
		totalHBondTypes = result[3]
		
		for molId in totalHBonds:
			
			avgHBonds = totalHBonds[molId] / framesCount
			
			outputFile.write("Molecule with id {} has an average of {} hydrogen bonds across the simulation.\n".format(molId, avgHBonds))
		
		outputFile.write("\n---------- HBonds Chains ----------\n")
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
		
		outputFile.write("\n---------- Total Number of HBonds ----------\n")
		
		for hbondType in totalHBondTypes:
			outputFile.write("HBond Type {} was found {} time(s).\n".format(hbondType, totalHBondTypes[hbondType]))
		
		outputFile.flush()
	

def main():
	
	try:
		# sys.argv[0] = script name.
		inputFileName = sys.argv[1]
	except:
		inputFileName = "res/input.txt"
	
	atomsFileName, atomsPerMolecule, fromFrame, toFrame, maxAngle, maxHBondDistance, maxBondDistance, hbondTypes = loader.loadInput(inputFileName)
	
	startTime = time.time()
	
	frames = loader.loadAnimatedXyz(atomsFileName, atomsPerMolecule, fromFrame=fromFrame, toFrame=toFrame)
	
	timeToLoadFrames = time.time() - startTime
	
	minutesToLoadFrames = timeToLoadFrames // 60
	secondsToLoadFrames = timeToLoadFrames % 60
	
	framesCount = len(frames)
	
	print("Loaded {} frames in {minuteMessage}{} seconds.".format(framesCount, secondsToLoadFrames, minuteMessage=(str(minutesToLoadFrames) + " minutes and") if minutesToLoadFrames > 0 else ""))
	
	startTime = time.time()
	
	result = analyzeFrames(frames, np.radians(maxAngle), maxHBondDistance, maxBondDistance, hbondTypes)
	
	timeToCountHBonds = time.time() - startTime
	
	try:
		outputFileName = sys.argv[2]
	except:
		outputFileName = "res/output.txt"
	
	outputResult(result, framesCount, outputFileName, maxHBondDistance, maxAngle)
	
	minutesToCountHBonds = timeToCountHBonds // 60
	secondsToCountHBonds = timeToCountHBonds % 60
	
	print("Took {minuteMessage}{} seconds to analyze {} frames; that's {} frames/second.".format(secondsToCountHBonds, framesCount, framesCount / timeToCountHBonds, minuteMessage=(str(minutesToCountHBonds) + " minutes and ") if minutesToCountHBonds > 0 else ""))

if __name__ == "__main__":
	main()
