from msmbuilder.hmm import GaussianHMM as GHMM
import numpy as np
import math
from scipy.spatial.distance import mahalanobis as mahadist
from copy import deepcopy



def hellinger(mean1, var1, mean2, var2):
#
	det1 = 1.0
	det2 = 1.0
	det12 = 1.0
	
	for i in range(len(var1)):
	#
		det1 *= var1[i]
		
		det2 *= var2[i]
		
		det12 *= 0.5*(var1[i] + var2[i])
	#
	
	term1 = math.sqrt(det1)*math.sqrt(det2)/det12
	
	term2 = 0.0
	
	for i in range(len(var1)):
	#
		term2 += -0.125*(mean1[i] - mean2[i])**2*2.0/(var1[i] + var2[i])
	#
	
	return math.sqrt(1-math.sqrt(term1)*math.exp(term2))
#



def prob_element(frame, means, ivars, popls, cluster):
#
	dividend = popls[cluster] * np.exp(0.5 * mahadist(frame, means[cluster], ivars[cluster])**2) / math.sqrt(((2*math.pi)**len(means[cluster])) * np.linalg.det(ivars[cluster]))
	
	divisor = 0.0
	
	for i in range(len(popls)):
	#
		divisor += popls[i] * np.exp(0.5 * mahadist(frame, means[i], ivars[i])**2) / math.sqrt(((2*math.pi)**len(means[i])) * np.linalg.det(ivars[i]))
	#
	
	return dividend / divisor
#


########################################################################
#																	   #
#	 		  					Node class							   #
#																	   #
########################################################################


class Node:
#
	def __init__(self, parent, location):
	#
		self.data = [] # .. # The leaves
		self.next = [] # .. # The child nodes
		self.parent = parent # .. # Location of the parent node (single inheritance)
		self.location = location # .. # Location of the node in the tree
	#
	
	def add_next(self, next):
	#
		#print("Adding next node...")
		
		location = deepcopy(self.location)
		
		location.append(len(self.next))
		
		#print("Parent location :")
		#print(self.location)
		#print("Node location :")
		#print(location)
		
		if next is None:
		#
			new_node = Node(self.location, location)  # .. # Creates a new node to append
			
			self.next.append(new_node)
		#
		else:
		#
			next.parent = self.location
			next.location = location
			
			self.next.append(next) # .. # Appends an existing node
		#
	#
	
	def add_data(self, data): # .. # Appends a single to data to existing data
	#
		self.data.append(data)
	#
	
	def extend_data(self, data): # .. # Appends a list of data to existing data
	#
		self.data.extend(data)
	#
	
	def print_node(self):
	#
		level = ""
		
		for i in range(len(self.location)):
		#
			level += "\t|"
		#
		
		print(level)
		
		if not (self.parent is None):
		#
			print("{:s}--> parent : {:s}".format(level, ", ".join(map(str, self.parent))))
		#
		else:
		#
			print("{:s}--> Head".format(level))
		#
		print("{:s}--> location : {:s}".format(level, ", ".join(map(str, self.location))))
		
		print("{:s}--> data : {:s}".format(level, ", ".join(map(str, self.data))))
		
		for i in range(len(self.next)):
		#
			self.next[i].print_node()
		#
	#
	
	def move(self):
	#
		location = deepcopy(self.location)
		
		location.append(0)
		
		for i in range(len(self.next)):
		#
			location[-1] = i
			
			self.next[i].location = location
			
			self.next[i].parent = self.location
			
			self.next[i].move()
		#
	#
	
	def search_data(self, data, found):
	#
		for i in self.data:
		#
			if i == data:
			#
				found.append(self.location)
			#
		#
		
		for node in self.next:
		#
			node.search_data(data, found)
		#
	#
	
	def remove(self): # .. # Destructor of node class
	#
		while self.next: # .. # Resursive deletion of all child nodes
		#
			self.next[-1].remove()
			
			del self.next[-1]
		#
		
		self.data = []
		self.parent = []
		self.location = []
		
		del self.next
		del self.data
		del self.parent
		del self.location
	#
#


########################################################################
#																	   #
#	 		  					Tree class							   #
#																	   #
########################################################################


class Tree:
#
	def __init__(self):
	#
		self.head = Node(None, [])
	#
	
	def get_node(self, location):
	#
		if not location:
		#
			return self.head
		#
		else:
		#
			current = self.head.next[location[0]]
			
			for i in range(1,len(location)):
			#
				current = current.next[location[i]]
			#
			
			return current
		#
	#
	
	def deepen(self, location, levels):
	#
		current = self.get_node(location)
		
		for i in range(levels):
		#
			current.add_next(None)
			
			current = current.next[-1]
		#
	#
	
	def broaden(self, location, n_new):
	#
		current = self.get_node(location)
		
		for i in range(n_new):
		#
			current.add_next(None)
		#
	#
	
	def add_data(self, location, data):
	#
		current = self.get_node(location)
		
		current.add_data(data)
	#
	
	def move_node(self, location1, location2):
	#
		cargo = self.get_node(location1)
		
		dest = self.get_node(location2)
		
		parent = self.get_node(cargo.parent)
		
		dest.add_next(cargo)
		
		cargo.move()
		
		del parent.next[location1[len(location1) - 1]]
		
		parent.move()
	#
	
	def move_data(self, location1, pos, location2):
	#
		giver = self.get_node(location1)
		
		taker = self.get_node(location2)
		
		data = giver.data[pos]
		
		del giver.data[pos]
		
		taker.add_data(data)
	#
	
	def move_all_data(self, location1, location2):
	#
		giver = self.get_node(location1)
		
		taker = self.get_node(location2)
		
		data = deepcopy(giver.data)
		
		taker.extend_data(data)
		
		del data
		
		giver.data = []
	#
	
	def delete_node(self, location):
	#
		current = self.get_node(location)
		
		parent = self.get_node(current.parent)
		
		pos = deepcopy(current.location[-1])
		
		current.remove()
		
		del parent.next[pos]
		
		for i in range(pos, len(parent.next)):
		#
			parent.next[i].location[-1] -= 1
			
			parent.next[i].move()
		#
	#
	
	def delete_data(self, location, pos):
	#
		current = self.get_node(location)
		
		del current.data[pos]
	#
	
	def wipe_data(self, location):
	#
		current = self.get_node(location)
		
		current.data = []
	#
	
	def print_tree(self):
	#
		self.head.print_node()
	#
	
	def search_data(self, data):
	#
		if not (isinstance(data, int), isinstance(data, float), isinstance(data, str)):
		#
			raise Exception("Cannot search for given type of data")
		#
		else:
		#
			found = []
			
			self.head.search_data(data, found)
			
			return found
		#
	#
#


########################################################################
#																	   #
#	 		  				Inheritance class						   #
#																	   #
########################################################################


class inheritance:
#
	def __init__(self, means, sigmas, popls, transmat):
	#
		self.means = means
		self.sigmas = sigmas
		self.pops = popls
		self.trans = transmat
		self.superstates = Tree()
		self.node = dict()
		self.lumps = []
		
		self.superstates.deepen(None, len(means[0]) - 1)
		
		loc = [0 for i in range(len(means[0]) - 1)]
		
		self.superstates.broaden(loc, len(means))
		
		loc.append(-1)
				
		for i in range(len(means)):
		#
			loc[-1] = i
			
			self.superstates.add_data(loc, i)
			
			self.node[i] = self.superstates.get_node(loc)
		#
		
		self.superstates.print_tree()
	#
	
	def same_lump(self, state1, state2):
	#
		if len(self.node[state1].location) == len(self.node[state2].location):
		#
			for i in range(len(self.node[state1].location)):
			#
				if self.node[state1].location[i] != self.node[state2].location[i]:
				#
					return False
				#
			#
		#
		else:
		#
			return False
		#
		
		return True
	#
	
	def lump(self):
	#
		for i in range(len(self.means)-1):
		#
			for j in range(i+1, len(self.means)):
			#
				if not self.same_lump(i,j):
				#
					prob = hellinger(self.means[i], self.sigmas[i], self.means[j], self.sigmas[j])
					
					#print("Hellinger distance between state {:d} and state {:d} : {:f}".format(i, j, prob))
					
					if prob < 0.95:
					#
						#print("Lumping states {:d} and {:d}...".format(i,j))
						
						self.superstates.move_all_data(self.node[j].location, self.node[i].location)
						
						to_del = deepcopy(self.node[j].location)
						
						for k in self.node[i].data:
						#
							self.node[k] = self.node[i]
						#
						
						self.superstates.delete_node(to_del)
						
						del to_del
					#
				#
			#
		#
		
		print("Lumping over :\n")
		
		self.superstates.print_tree()
		
		print("\n")
	#
	
	def match_lumps(self, node1, node2, means, sigmas):
	#
		match = False
		
		for i in node1.data:
		#
			for j in node2.data:
			#
				prob = hellinger(means[i], sigmas[i], means[j], sigmas[j])
				
				if prob < 0.95:
				#
					match = True
				#
			#
		#
		
		if not match:
		#
			for node_1 in node1.next:
			#
				for node_2 in node2.next:
				#
					match = self.match_lumps(node_1, node_2, means, sigmas)
				#
			#
		#
				
		if not match:
		#
			for node in node2.next:
			#
				match = self.match_lumps(node1, node, means, sigmas)
			#
		#
		
		if not match:
		#
			for node in node1.next:
			#
				match = self.match_lumps(node, node2, means, sigmas)
			#
		#
		
		return match
	#
	
	def retrolump(self):
	#
		means = self.means
		sigmas = self.sigmas
		
		parent_loc = self.node[0].location
		
		while len(means[0]) > 1:
		#
			for mean in means:
			#
				del mean[-1]
			#
			
			for sigma in sigmas:
			#
				del sigma[-1]
			#
			
			del parent_loc[-1]
			
			parent_node = self.superstates.get_node(parent_loc)
			
			grandparent_node = self.superstates.get_node(parent_node.parent)
			
			while parent_node.next:
			#
				self.superstates.broaden(grandparent_node.location, 1)
				
				new_parent = grandparent_node.next[-1]
				
				self.superstates.move_node(parent_node.next[0].location, new_parent.location)
				
				i = -1
				
				while i < len(new_parent.next) - 1:
				#
					i += 1
					
					for node in parent_node.next:
					#
						if self.match_lumps(new_parent.next[i], node, means, sigmas):
						#
							self.superstates.move_node(node.location, new_parent.location)
						#
					#
				#
			#
			
			self.superstates.delete_node(parent_loc)
			
			print("\n\nRetrolumping done on this layer :\n")
			
			self.superstates.print_tree()
			
			print("\n")
		#
		
		self.superstates.print_tree()
	#
#
