'''
Created on Jun 30, 2017

@author: Rawad
'''

class Graph:
	'''
	Uses an adjacency list representation to represent an undirected graph.
	'''
	
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
	
	def __len__(self):
		return len(self.graph)
	
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
	
	def getConnectedComponents(self):
		'''
		Returns list of graphs, each representing a connected component cotnained in this original graph.
		'''
		
		visited = {i: False for i in self.graph}
		
		graphs = []
		
		for v in self.graph:
			if not visited[v]:
				graph = Graph()
				graphs.append(graph)
				self.connectedComponentsUtil(v, visited, graph)
		
		return graphs
	
	def connectedComponentsUtil(self, v, visited, graph):
		
		visited[v] = True
		
		for u in self.graph[v]:
			if not visited[u]:
				graph.addGroup(v, u)
				self.connectedComponentsUtil(u, visited, graph)
