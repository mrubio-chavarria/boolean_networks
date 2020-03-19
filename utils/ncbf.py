
"""
INDICATIONS:
	In this file, they are given several utilities to the calculation of nested canalizing
boolean functions (boolean_networks).
"""
from itertools import combinations
from functools import reduce
from operator import add
from random import choice
from numpy import linspace


def get_level(tree, path, level):
	"""
	DESCRIPTION:
	:param tree: list of lists with all possible groups to be formed.
	:param path: list with the groups so far.
	:param level: the level in the tree from we are to get the available groups.
	:return: layers: the list with all the available groups.
	"""
	packed_path = ''.join(path)
	candidates = tree[level]
	layers = [candidate for candidate in candidates if not any([True if letter in packed_path else False for letter
																in candidate])]
	return layers

def charCal(row):
	"""
	DESCRIPTION:
	:param row: list of string with different length.
	:return: weight: number of letters in the whole list.
	"""
	weight = sum([len(element) for element in row])
	return weight

def pathCalc(row, paths):
	"""
	DESCRIPTION: 
	A function to calculate a path in the tree of functions
	"""

	# Parameters
	act = row.loc['activators'] or []
	inh = row.loc['inhibitors'] or []
	e_act = len(act)
	e_inh = len(inh)
	p_total = e_act + e_inh

	# Tree of activators
	tree_act = [act]
	for i in range(2, len(act)+1):
		level = tree_act.append([reduce(add, tup) for tup in list(combinations(act, i)) if tup is not None])
		if level is not None:
			tree_act.append(level)

	# Tree of inhibitors
	tree_inh = [inh]
	for i in range(2, len(inh)+1):
		level = tree_inh.append([reduce(add, tup) for tup in list(combinations(inh, i)) if tup is not None])
		if level is not None:
			tree_inh.append(level)

	# Start to calculate the paths
	paths = []
	# Develop a set of paths
	flat_tree = [item for sublist in tree_act for item in sublist]
	for initial_leaf in flat_tree:
		# Initialize the parameters
		path = [initial_leaf]
		var_act = e_act - len(initial_leaf)
		var_inh = e_inh
		p_ac = len(initial_leaf)
		p_c = p_total - p_ac - var_act
		while True:
			# Set the level in the tree
			if var_act == 0:
				next_level = p_c - 1
			else:
				next_level = int(choice(linspace(0, p_c - 1, p_c)))
			# Add inhibitor layer
			local_inh = get_level(tree_inh, path, next_level)
			next_leaf = choice(local_inh)
			path.append(next_leaf)
			# Update parameters
			p_ac = charCal(path)
			var_inh += -len(next_leaf)
			p_c = p_total - p_ac - var_inh
			# Assess the exit
			if p_ac == p_total:
				break

			# Set the level in the tree
			if var_inh == 0:
				level = p_c - 1
			else:
				level = int(choice(linspace(0, p_c - 1, p_c)))
			# Add activator layer
			local_act = get_level(tree_act, path, level)
			leaf = choice(local_act)
			path.append(leaf)
			# Update parameters
			p_ac = charCal(path)
			var_act += -len(leaf)
			p_c = p_total - p_ac - var_act
			# Assess the exit
			if p_ac == p_total:
				break

		paths.append(path)
		rev = path[::-1]
		paths.append(rev)
	return paths

def ncbfCalc(data):
	"""
	DESCRIPTION: 
	With this function, given a graph, we obtain all the possible NCBFs. 
	"""
	paths = []
	for index, row in data.iterrows():
		# Obtain the number of possible layers
		paths.append(pathCalc(row=row, paths=[]))
	return paths
