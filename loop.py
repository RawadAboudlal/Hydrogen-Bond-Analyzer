'''
Created on Jun 30, 2017

@author: Rawad
'''

class Graph:
	
	def __init__(self):
		self.graph = {}
	
	def addGroup(self, a, b):
		
		if not a in self.graph:
			self.graph[a] = []
		
		if not b in self.graph:
			self.graph[b] = []
		
		if not a in self.graph[b]:
			self.graph[b].append(a)
		
		if not b in self.graph[a]:
			self.graph[a].append(b)
	
	def isCyclic(self):
		
		visited = {i: False for i in self.graph}
		
		for v in self.graph:
			if not visited[v]:
				if self.isCyclicUtil(visited, v, -1):
					return True
				
		
		return False
	
	def isCyclicUtil(self, visited, v, parent):
		
		visited[v] = True
		
		for u in self.graph[v]:
			if not visited[u]:
				if self.isCyclicUtil(visited, u, v):
					return True
			elif parent != u:
				return True
				
		return False

def createGraphFromHBonds(hbondList):
	
	graph = Graph()
	
	for hbond in hbondList:
		graph.addGroup(hbond.mol1, hbond.mol2)
	
	return graph
