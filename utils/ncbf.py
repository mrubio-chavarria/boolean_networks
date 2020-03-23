
"""
INDICATIONS:
	In this file, they are given several utilities to the calculation of nested canalizing
boolean functions (boolean_networks).
"""
from itertools import combinations, permutations, product, groupby
from functools import reduce
from operator import add
from random import choice
from numpy import linspace
import numpy as np


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
    :return: number of letters in the whole list.
    """
    return sum([len(element) for element in row])

def relatedPaths(path, path_params, tree, other_tree):
    """
    DESCRIPTION:
    A function that, given a path, returns all related paths.
    :param path: a set of layers.
    :param path_params: calculations of each step of the path.
    :param tree: tree of the inital step of the path.
    :return: set of all related paths.
    """
    # Parameters
    params = np.array(path_params)
    paths = []
    # All combinations with the same structure
    even_perms = list(permutations(path[::2]))
    odd_perms = list(permutations(path[1::2]))
    perms = list(product(even_perms, odd_perms))
    for perm in perms:
        path[::2] = list(perm[0])
        path[1::2] = list(perm[1])
        paths.append(path)
    return paths

def auxiliar_function(groupA, groupB, n, path=None):
    """
    DESCRIPTION:
    A generator to iterate over a dynamic number of for loops. Here are generated all possible structures given a set of
    groups.
    :param groupA: first set of levels in the structure.
    :param groupB: second set of levels in the structure.
    :param n: total weight of the path.
    :param path: path to propagate, to generate.
    :return: all possible structures, paths.
    """
    # Parameters to assess the path
    n_A = groupA[-1]
    n_B = groupB[-1]
    if path:
        if len(path)%2 == 0:
            # Last layer A, add B
            for el_B2 in groupB:
                pathB2 = path + [el_B2]
                valB = np.sum(pathB2)
                if valB == n:
                    model = np.array(pathB2)
                    vecA = np.zeros([1, len(model)])
                    vecA[0, 0::2] = 1
                    vecB = np.zeros([1, len(model)])
                    vecB[0, 1::2] = 1
                    if (np.sum(vecB * model) == n_B) and (np.sum(vecA * model) == n_A):
                        yield pathB2
                elif valB < n:
                    yield from auxiliar_function(groupA, groupB, n, path=pathB2)
        else:
            # Last layer B, add A
            for el_A2 in groupA:
                pathA2 = path + [el_A2]
                valA = np.sum(pathA2)
                if valA == n:
                    model = np.array(pathA2)
                    vecA = np.zeros([1, len(model)])
                    vecA[0, 0::2] = 1
                    vecB = np.zeros([1, len(model)])
                    vecB[0, 1::2] = 1
                    if (np.sum(vecB * model) == n_B) and (np.sum(vecA * model) == n_A):
                        yield pathA2
                elif valA < n:
                    yield from auxiliar_function(groupA, groupB, n, path=pathA2)
    else:
        # Start with a layer of A
        for el_A1 in groupA:
            pathA1 = [el_A1]
            yield from auxiliar_function(groupA, groupB, n, path=pathA1)

def structsCalc(groupA, groupB, n):
    """
    DESCRIPTION:
    A function that achieves all the combinations of numbers from 0 to groupA, and from 0 to groupB in an iterative way
    and each row sums n.
    :param groupA: number of elements in the first group.
    :param groupB: number of elements in the second group.
    :param n: total value by group.
    :return: list of possible structures.
    """
    groupA = np.linspace(1, groupA, groupA)
    groupB = np.linspace(1, groupB, groupB)
    selectedA = list(auxiliar_function(groupA, groupB, n))
    selectedB = list(auxiliar_function(groupB, groupA, n))
    return [selectedA, selectedB]

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

    # Calculate combinations of levels
    structs = structsCalc(len(tree_act), len(tree_inh), p_total)

    # Start to calculate the paths
    paths = []
    paths_params = {}
    count = 0
    for branch in tree_act:
        # Initialize the parameters
        initial_leaf = branch[0]
        path = [initial_leaf]
        path_params = []
        var_act = e_act - len(initial_leaf)
        var_inh = e_inh
        p_ac = len(initial_leaf)
        p_c = p_total - p_ac - var_act
        level = 0
        # Store the step parameters
        path_params.append([p_total, p_ac, var_act, var_inh, p_c, level])
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
            # Store the step parameters
            path_params.append([p_total, p_ac, var_act, var_inh, p_c, next_level])
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
            # Store the step parameters
            path_params.append([p_total, p_ac, var_act, var_inh, p_c, level])
            # Assess the exit
            if p_ac == p_total:
                break
        # Set of related paths
        paths_set = relatedPaths(path, path_params, tree_act, tree_inh)
        # Store path and param
        paths = paths + paths_set
        paths_params[count] = path_params
        count += 1

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
