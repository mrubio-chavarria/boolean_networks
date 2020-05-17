
"""
INDICATIONS:
In this file, they are given several utilities to the calculation of nested canalizing boolean functions
(boolean_networks).
"""

from itertools import combinations, permutations, product
from functools import reduce
from operator import add
from random import choice
import pandas as pd
import numpy as np
from utils.stp import lGen
from utils.stp import stpn, stp


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
    new_path = np.zeros([1, len(path)]).tolist()[0]
    params = np.array(path_params)
    paths = []
    # All combinations with the same structure
    even_perms = list(permutations(path[::2]))
    odd_perms = list(permutations(path[1::2]))
    perms = list(product(even_perms, odd_perms))
    for perm in perms:
        new_path[::2] = list(perm[0])
        new_path[1::2] = list(perm[1])
        if new_path in paths:
            continue
        paths.append(new_path[:])
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

def pathCalc(row):
    """
    DESCRIPTION:
    A function to calculate a path in the tree of functions
    :param row: a row of the networks table. Activators and inhibitors.
    :return: a list with all possible functions which meet with the row element.
    """

    # Parameters
    act = row.loc['activators'] or []
    inh = row.loc['inhibitors'] or []
    e_act = len(act)
    e_inh = len(inh)
    p_total = e_act + e_inh

    # Special cases of empty fields
    if act == [''] and inh != ['']:
        return [[''.join(inh)]]
    elif act != [''] and inh == ['']:
        return [[''.join(act)]]
    elif act == [''] and inh == ['']:
        return ['INPUT']

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

    # Calculate the paths
    paths = []
    paths_params = {}
    for s in range(0, 2):
        if s == 0:
            even_tree = tree_act
            e_even = e_act
            odd_tree = tree_inh
            e_odd = e_inh
        else:
            even_tree = tree_inh
            e_even = e_inh
            odd_tree = tree_act
            e_odd = e_even
        structs_set = structs[s]

        count = 0
        for struct in structs_set:

            # Initialize the parameters
            struct = np.array(struct) - 1
            struct = struct.astype(int).tolist()
            initial_leaf = choice(even_tree[struct[0]])
            path = [initial_leaf]
            path_params = []
            var_even = e_even - len(initial_leaf)
            var_odd = e_odd
            p_ac = len(initial_leaf)
            p_c = p_total - p_ac - var_even
            level = 0

            # Store the step parameters
            path_params.append([p_total, p_ac, var_even, var_odd, p_c, level])
            struct_count = 1
            while True:
                # Set the level in the tree
                if var_even == 0:
                    next_level = p_c - 1
                else:
                    # Old by chance version:
                    # next_level = int(choice(linspace(0, p_c - 1, p_c)))
                    next_level = struct[struct_count]
                    struct_count += 1
                # Add inhibitor layer
                local_odd = get_level(odd_tree, path, next_level)
                next_leaf = choice(local_odd)
                path.append(next_leaf)
                # Update parameters
                p_ac = charCal(path)
                var_odd += -len(next_leaf)
                p_c = p_total - p_ac - var_odd
                # Store the step parameters
                path_params.append([p_total, p_ac, var_even, var_odd, p_c, next_level])
                # Assess the exit
                if p_ac == p_total:
                    break

                # Set the level in the tree
                if var_even == 0:
                    level = p_c - 1
                else:
                    # Old by chance version:
                    # level = int(choice(linspace(0, p_c - 1, p_c)))
                    level = struct[struct_count]
                    struct_count += 1
                # Add activator layer
                local_even = get_level(even_tree, path, level)
                leaf = choice(local_even)
                path.append(leaf)
                # Update parameters
                p_ac = charCal(path)
                var_even += -len(leaf)
                p_c = p_total - p_ac - var_even
                # Store the step parameters
                path_params.append([p_total, p_ac, var_even, var_odd, p_c, level])
                # Assess the exit
                if p_ac == p_total:
                    break
            # Set of related paths
            paths_set = relatedPaths(path, path_params, even_tree, odd_tree)
            # Store path and param
            paths = paths[:] + paths_set[:]
            paths_params[count] = path_params
            count += 1

    return paths

def ncbfCalc(data):
    """
    DESCRIPTION:
    With this function, given a graph, we obtain all the possible NCBFs.
    :param data: a dataframe with the representation of the network. One row per node, first column for activators and
    second for inhibitors.
    :return: all possible NCBF for each gene according with the structure presented in Murrugarra 2013 for NCBF.
    """
    paths = []
    for index, row in data.iterrows():
        # Obtain the number of possible layers
        paths.append(pathCalc(row=row))

    # Build the response dataframe
    paths = pd.Series(data=paths, index=data.index)
    return paths

def networksCalc(paths, path=None, index=0):
    """
    DESCRIPTION:
    Given an object of all possible combinations, the one offered by the function above, we build another with all
    possible networks. One network per line, one node per column.
    :param data: the set of all possible combinations of functions.
    :return: the network. Due to the fact that we have a generator, it will be a list of objects.
    """
    if path is not None:
        if index%2 == 0:
            for step in paths.iloc[index]:
                path_second = path[:] + [step]
                index_second_A = index + 1
                if index_second_A >= len(paths.index):
                    yield path_second
                else:
                    yield from networksCalc(paths, path=path_second, index=index_second_A)
        else:
            for step in paths.iloc[index]:
                path_second = path[:] + [step]
                index_second_B = index + 1
                if index_second_B >= len(paths.index):
                    yield path_second
                else:
                    yield from networksCalc(paths, path=path_second, index=index_second_B)
    else:
        index = 0
        for step in paths.iloc[index]:
            path_first = [step]
            index_first = index + 1
            yield from networksCalc(paths, path=path_first, index=index_first)


def netValidator(networks, graph, attractors):
    """
    DESCRIPTION:
    Given a function and an attractor, it is returned whether the function meets the attractors condition or not.
    :param networks: [list] NCBFs to be validated with the structure of Murrugarra 2013.
    :param graph: [pandas dataframe] graph from which come the networks in the structure set in previous stages.
    :return: [list] networks that have passed the validation.
    """
    # Prepare the attractors
    t = np.array([[1], [0]])
    f = np.array([[0], [1]])
    for i in range(0, len(attractors[:])):
        attractors[i] = stpn([t if num == 1 else f for num in attractors[i]])

    # Execute the validation
    final_networks = []
    for net in networks:
        l = lGen(net, graph)
        validation = all([True if all(attractor == stp(l, attractor)) else False for attractor in attractors])
        if validation:
            final_networks.append(net)

    return final_networks


