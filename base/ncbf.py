#! /usr/bin/env python3
# -*- coding: utf-8 -*-

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
from alive_progress import alive_bar
from base.exceptions import InputAlterationException, ConvergenceException
from PyBoolNet import StateTransitionGraphs, Attractors, FileExchange
import itertools
import random
from graphs.conflicts_theory import Pathway, KMap, ConflictsManager
from base.hibrids import Result
import progressbar


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
    if path is not None:
        if len(path) % 2 != 0:
            # Last layer A, add B
            for el_B2 in groupB[:]:
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
            for el_A2 in groupA[:]:
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
        for el_A1 in groupA[:]:
            yield from auxiliar_function(groupA, groupB, n, path=[el_A1])

def structsCalc(groupA, groupB, n):
    """
    DESCRIPTION:
    A function that achieves all the combinations of numbers from 0 to groupA, and from 0 to groupB in an iterative way
    and each row sums n.
    :param groupA: [int] number of elements in the first group.
    :param groupB: [int] number of elements in the second group.
    :param n: [int] total value by group.
    :return: [list] all possible structures.
    """
    groupA = np.linspace(1, groupA, groupA)
    groupB = np.linspace(1, groupB, groupB)
    selectedA = list(auxiliar_function(groupA, groupB, n))
    selectedB = list(auxiliar_function(groupB, groupA, n))
    return [selectedA, selectedB]

def pathCalc(row, tags):
    """
    DESCRIPTION:
    A function to calculate a path in the tree of functions
    :param row: a row of the networks table. Group 1 and 2 (Activators and inhibitors).
    :return: a list with all possible functions which meet with the row element.
    """

    # Parameters
    act = row.loc[tags[0]] or []
    inh = row.loc[tags[1]] or []
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


def ncbfCalc(data, tags):
    """
    DESCRIPTION:
    With this function, given a graph, we obtain all the possible NCBFs.
    :param data: [pandas DataFrame] with the representation of the network. One row per node, first column for activators and
    second for inhibitors.
    :return: [pandas Series] all possible NCBF for each gene according with the structure presented in Murrugarra 2013 for NCBF.
    """
    paths = []
    for index, row in data.iterrows():
        # Obtain the number of possible layers
        paths.append(pathCalc(row=row, tags=tags))

    # Build the response dataframe
    paths = pd.Series(data=paths, index=data.index)
    return paths

def conflict_ncbfCalc(variables, tag):
    """
    DESCRIPTION:
    It does the same as the function above but with all possible combinations in groups of the given variables.
    :param variables: [pandas DataFrame] variables of the network. One row per node.
    :return: [tuple] all possible NCBF for each gene according with the structure presented in Murrugarra 2013 for NCBF.
    """
    # Auxiliary functions
    def aux_i(node):
        params = variables.loc[node, tag]
        for i in range(0, len(params) + 1):
            yield itertools.combinations(params, i)

    # Parameters
    nodes = variables.index

    # Base dataframe
    index = nodes
    columns = ['activators', 'inhibitors']

    # Make all combinations
    all_combs = []
    for node in nodes:
        combs = [item for sublist in list(aux_i(node)) for item in sublist]
        combs = [comb for comb in list(itertools.product(combs, repeat=2))
                 if len(comb[0] + comb[1]) == len(variables.loc[node, tag])
                 if not any([True if element in comb[1] else False for element in comb[0]])]
        combs = [[group if len(group) > 0 else [''] for group in comb] for comb in combs]
        all_combs.append(combs)
    all_combs = itertools.product(*all_combs)
    base = lambda x: pd.DataFrame(data=x, index=index, columns=columns)
    frames = list(map(base, all_combs))

    # Build the paths
    all_paths = []
    for frame in frames:
        paths = ncbfCalc(frame, columns)
        all_paths.append(paths)

    return all_paths, frames

def networksCalc(paths, path=None, index=0):
    """
    DESCRIPTION:
    Given an object of all possible combinations, the one offered by the function above, we build another with all
    possible networks. One network per line, one node per column.
    :param data: the set of all possible combinations of functions.
    :return: the network. Due to the fact that we have a generator, it will be a list of objects.
    """
    if path is not None:
        if index % 2 == 0:
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

def conflict_net2boolnet(net, graph, tags, groups=None):
    """
    DESCRIPTION:
    """

    # Auxiliary functions
    def aux_I():
        for node in groups.keys():
            group1 = groups[node][list(groups[node].keys())[0]]
            group1 = group1 if len(group1) > 0 else [()]
            group2 = groups[node][list(groups[node].keys())[1]]
            group2 = group2 if len(group2) > 0 else [()]
            new_groups = list(itertools.product(group1, group2, repeat=1))
            for group in new_groups:
                expr = network[node]
                variables = [item for sublist in group for item in sublist]
                # Due to the properties of the NCBF we can assume that there will not be any variable twice
                positions = [expr.index(variable[0]) for variable in variables if variable[0] in expr]
                i = 0
                while i < len(expr):
                    if i in positions:
                        p = positions.index(i)
                        value = variables[p][1]
                        variables.pop(p)
                        positions.pop(p)
                        if value == 0:
                            expr = expr[0:i - 1] + expr[i::]
                            positions = [pos - 1 for pos in positions if pos > i]
                            i += -1
                    i += 1
                yield node, expr

    def aux_II(groups, counts):
        """
        DESCRIPTION:
        A function to perform an absolutely customized cartesian product.
        """
        selections = [combinations(g, c) for g, c in zip(groups, counts)]
        for n_tuple in product(*selections):
            yield tuple(itertools.chain.from_iterable(n_tuple))

    # Parameters
    activators = graph[tags[0]]
    tags = list(graph.index)
    network = {}

    # Section to manage the original nodes
    for i in range(0, len(net)):
        node = net[i]
        header = ''
        tail = ''
        if node == 'INPUT':
            network[tags[i]] = tags[i]
            continue
        if node[0][0] in activators[tags[i]]:
            header = '!('
            tail = ')'
        last = True
        load = ''
        for layer in reversed(node):
            if last:
                load = '&'.join(list(map(lambda x: '!' + x, layer)))
                last = False
                continue
            local_load = '&'.join(list(map(lambda x: '!' + x, layer)))
            load = f'{local_load}&!({load})'
        node = header + load + tail
        network[tags[i]] = node

    # Section to tackle the conflict nodes
    if groups is not None:
        combs = list(aux_I())
        groups = {key: [] for key in tags}
        list(map(lambda x: groups[x[0]].append(x[1]), combs))
        groups = [groups[key] for key in groups.keys()]
        counts = [1] * len(groups)
        combs = list(aux_II(groups, counts))
        for comb in combs:
            network = {tags[i]: comb[i] for i in range(0, len(tags))}
            yield network

def net2boolnet(net, graph, tags):
    """
    DESCRIPTION:
    A function to parse the incoming network to a Boolnet format in string. The network to be taken by PyBoolNet.
    :param net: [list] the network in our format.
    :return: [string] the network in Boolnet format.
    """

    # Parameters
    activators = graph[tags[0]]
    tags = list(graph.index)
    network = {}

    # Section to manage the original nodes
    for i in range(0, len(net)):
        node = net[i]
        header = ''
        tail = ''
        if node == 'INPUT':
            network[tags[i]] = tags[i]
            continue
        if node[0][0] in activators[tags[i]]:
            header = '!('
            tail = ')'
        last = True
        load = ''
        for layer in reversed(node):
            if last:
                load = '&'.join(list(map(lambda x: '!' + x, layer)))
                last = False
                continue
            local_load = '&'.join(list(map(lambda x: '!' + x, layer)))
            load = f'{local_load}&!({load})'
        node = header + load + tail
        network[tags[i]] = node

    network = '\n'.join([f'{item[0]}, {item[1]}' for item in network.items()])
    return network


def net2file(net, filename):
    """
    DESCRIPTION:
    A function to store the dict with the network in a file with boolnet format.
    :param net: [dict] network in boolnet format.
    :return: None.
    """
    file = open(filename, 'w')
    message = '\n'.join([f'{item[0]}, {item[1]}' for item in net.items()])
    file.write(message)
    file.close()
    return filename


def conflicts_manager(net, conflicts_networks, conflicts_graph, tags):
    """
    DESCRIPTION:
    A function to introduce the expressions of the conflicts in the complete network.
    """
    # Get all possible combinations
    groups = {node: {column: 1 for column in tags} for node in conflicts_graph.index}
    for node in conflicts_graph.index:
        for column in tags:
            local_groups = [item for sublist in [[(letter, 1), (letter, 0)] for letter in conflicts_graph[column][node]]
                            for item in sublist if item[0] != '']
            local_combs = list(itertools.combinations(local_groups, len(conflicts_graph[column][node])))
            local_combs = [list(comb) for comb in local_combs]
            condition = all(True if len(group) == 1 else False for group in local_combs)
            if not condition:
                groups[node][column] = [group for group in local_combs if group[0][0] != group[1][0] if
                                        len(set([item[0] for item in group])) == len(conflicts_graph[column][node])]
            else:
                groups[node][column] = local_combs


    # Calculate the subnetworks
    for network in conflicts_networks:
        networks = [neto for neto in list(conflict_net2boolnet(network, conflicts_graph, groups=groups, tags=['activators', 'inhibitors']))]
        """
        networks = [[{key: value[2:-1] for (key, value) in network.items()}, network] for network
                    in list(conflict_net2boolnet(network, conflicts_graph, groups=groups, tags=['activators', 'inhibitors']))]
        networks = [item for sublist in networks for item in sublist]
        """
        networks_set = [net + '\n' + '\n'.join([f'{item[0]}, {item[1]}' for item in network.items()]) for network in networks]
        yield networks_set


def netValidator(initial_networks=None, initial_graph=None, original_networks=None, original_graph=None, tags=None,
                 unfixed_sets_conflicts_networks=None, unfixed_conflicts_graphs=None, attractors=None):
    """
    DESCRIPTION:
    Given a function and an attractor, it is returned whether the function meets the attractors condition or not.
    :param initial_networks: [list] with structures of the networks. Lists containing the set of steps which form the
    path in the tree algorithm.
    :param initial_graph: [pandas DataFrame] the representation of the original network.
    :param tags: [list] labels for the columns with different roles for the positions in the columns of the data.
    Example: activators and inhibitors.
    :param original_networks: [list] structures devoted to the extension of the network. They incorporate the new nodes
    but they only represent the expressions of the original nodes.
    :param original_graph: [pandas DataFrame] the same concept of original networks but extrapolated to the graph. The
    expressions of the original nodes with the new nodes.
    :param attractors: [list] strings representing the steady states of the network which we are looking for.
    :param unfixed_sets_conflicts_networks: [list] sets of paths, with different structures for the networks according
    to each set of path. But these are partial structures, they only describe expressions for the new nodes.
    :param unfixed_conflicts_graphs: [list] the same as before but for the graphs. All the possible graphs describing
    only the new nodes. Different graphs, different set of paths with all possible combinations.
    :return: [list] networks that have passed the validation.
    """
    # Auxiliary functions
    def first_validation(unfixed_sets_conflicts_networks, unfixed_conflicts_graphs, limit=None):
        """
        DESCRIPTION:
        The function to extend the networks with the given conflicts sets.
        :param unfixed_sets_conflicts_networks: [list] sets of paths, with different structures for the networks
        according to each set of path. But these are partial structures, they only describe expressions for the new
        nodes.
        :param unfixed_conflicts_graphs: [list] the same as before but for the graphs. All the possible graphs
        describing only the new nodes. Different graphs, different set of paths with all possible combinations.
        :param limit: [int] number of networks to be extended.
        :return: [dictionary] final network with its attractors.
        """
        print('Extending networks')
        num_sets_conflicts_networks = len(unfixed_sets_conflicts_networks) if limit is None else limit
        assessed_networks = []
        with alive_bar(len(original_networks)*num_sets_conflicts_networks) as bar:
            for i in range(0, num_sets_conflicts_networks):
                conflicts_networks = unfixed_sets_conflicts_networks[i]
                conflicts_graph = unfixed_conflicts_graphs[i]
                # Set progress evaluation
                for original_net in original_networks:
                    # Finish the construction of all possible networks with conflicts networks
                    original_net = net2boolnet(original_net, original_graph, tags=['activators', 'inhibitors'])
                    nets = [item for sublist in list(conflicts_manager(original_net, conflicts_networks, conflicts_graph, tags))
                            for item in sublist]
                    # Execute the validation of all possible networks
                    for net in nets:
                        if net not in assessed_networks:
                            primes = FileExchange.bnet2primes(net)
                            stg = StateTransitionGraphs.primes2stg(primes, "synchronous")
                            steady, cyclic = Attractors.compute_attractors_tarjan(stg)
                            steadies = [
                                [att for att in attractors if steady_att.startswith(att)] for steady_att in steady
                            ]
                            steadies = list(set([item for sublist in steadies for item in sublist]))
                            condition = True if len(steadies) == len(attractors) else False
                            # Assess whether the result is to be sent or not
                            if condition:
                                assessed_networks.append(net)  # TO BE IMPROVED
                                yield {
                                    'network': net,
                                    'steady': steady,
                                    'cyclic': cyclic
                                }
                # Update progress
                bar()
        print('Networks extension completed')

    # Extend the original networks and execute the validation
    extended_ncbf_networks = list(first_validation(unfixed_conflicts_graphs=unfixed_conflicts_graphs,
                                                   unfixed_sets_conflicts_networks=unfixed_sets_conflicts_networks))
    return extended_ncbf_networks


def post_process(unfixed_sets_conflicts_networks, unfixed_conflicts_graphs, original_networks, original_graph,
                 attractors, limit=None):
    """
    DESCRIPTION:
    The function to extend the networks with the given conflicts sets.
    :param unfixed_sets_conflicts_networks: [list] sets of paths, with different structures for the networks
    according to each set of path. But these are partial structures, they only describe expressions for the new
    nodes.
    :param unfixed_conflicts_graphs: [list] the same as before but for the graphs. All the possible graphs
    describing only the new nodes. Different graphs, different set of paths with all possible combinations.
    :param limit: [int] number of networks to be extended.
    :param original_networks: [list] structures devoted to the extension of the network. They incorporate the new nodes
    but they only represent the expressions of the original nodes.
    :param original_graph: [pandas DataFrame] the same concept of original networks but extrapolated to the graph. The
    expressions of the original nodes with the new nodes.
    :param attractors: [list] strings representing the steady states of the network which we are looking for.
    :return: [dictionary] final network with its attractors.
    """
    print('Extending networks')
    num_sets_conflicts_networks = len(unfixed_sets_conflicts_networks) if limit is None else limit
    assessed_networks = []
    with alive_bar(len(original_networks)*num_sets_conflicts_networks) as bar:
        for i in range(0, num_sets_conflicts_networks):
            conflicts_networks = unfixed_sets_conflicts_networks[i]
            conflicts_graph = unfixed_conflicts_graphs[i]
            # Set progress evaluation
            for original_net in original_networks:
                # Finish the construction of all possible networks with conflicts networks
                original_net = net2boolnet(original_net, original_graph, tags=['activators', 'inhibitors'])
                nets = [item for sublist in list(conflicts_manager(original_net, conflicts_networks, conflicts_graph, tags))
                        for item in sublist]
                # Execute the validation of all possible networks
                for net in nets:
                    if net not in assessed_networks:
                        primes = FileExchange.bnet2primes(net)
                        stg = StateTransitionGraphs.primes2stg(primes, "synchronous")
                        steady, cyclic = Attractors.compute_attractors_tarjan(stg)
                        steadies = [
                            [att for att in attractors if steady_att.startswith(att)] for steady_att in steady
                        ]
                        steadies = list(set([item for sublist in steadies for item in sublist]))
                        condition = True if len(steadies) == len(attractors) else False
                        # Assess whether the result is to be sent or not
                        if condition:
                            assessed_networks.append(net)  # TO BE IMPROVED
                            yield {
                                'network': net,
                                'steady': steady,
                                'cyclic': cyclic
                            }
            # Update progress
            bar()
    print('Networks extension completed')


def third_validation(graph, initial_networks, attractors, simulations=20, max_iterations=2000):
    """
    DESCRIPTION:
    The third mechanism ot treating the conflicts.
    :param graph: [pandas DataFrame] the original representation of the network.
    :param initial_networks: [list] with structures of the networks. Lists containing the set of steps which form the
    path in the tree algorithm.
    :param simulations: [int] number of analysis of the same network that are to be performed.
    :param max_iterations: [int] maximum number of iterations that the recursive algorithm is allowed to perform. The
    objective is to avoid non-convergent solutions.
    :param attractors: [list] strings representing the steady states of the network which we are looking for.
    :return: [list] results which the attractors condition.
    """
    # Auxiliary functions
    def str_gen(n):
        for i in range(0, n):
            yield '0'

    # All possible variables combinations
    print('Initializing parameters and pathways')
    combinations = [bin(i).split('b')[1] if len(bin(i).split('b')[1]) == len(graph.index) else
                    ''.join(str_gen(len(graph.index) - len(bin(i).split('b')[1]))) + bin(i).split('b')[1]
                    for i in range(0, 2 ** len(graph.index))]
    space = [{graph.index[i]: int(c[i]) for i in range(0, len(c))} for c in combinations]

    # Generate the initial set of pathways
    initial_pathways = []
    with progressbar.ProgressBar(max_value=len(graph.index)) as bar:
        for i in range(0, len(graph.index)):
            node = graph.index[i]
            [initial_pathways.append(Pathway(antecedent=act, consequent=node, activator=True, space=space))
             for act in graph['activators'][node] if act != '']
            [initial_pathways.append(Pathway(antecedent=inh, consequent=node, activator=False, space=space))
             for inh in graph['inhibitors'][node] if inh != '']
            bar.update(i)

    # Set initial parameters of the similations
    results = {'accepted': [], 'dismissed': []}

    # Get all inputs
    inputs = [el for el in graph.index if graph.loc[el, 'activators'] == [''] and
              graph.loc[el, 'inhibitors'] == ['']]

    # Launch the simulations
    print('Launching simulations')
    with progressbar.ProgressBar(max_value=simulations) as bar:
        for m in range(0, simulations):
            # Generate the matrix with the priorities for the simulation
            blosum = blosum_generator(graph)
            net_counter = -1
            for initial_network in initial_networks:
                net_counter += 1
                condition_result = True  # set the condition to append the result of the simulation
                pathways = initial_pathways[:]
                n_pathways = len(pathways)  # number of pathways
                # Generate the basis map
                base_map = KMap(initial_network, graph)
                # Control parameters
                i = 0
                iter_count = 0
                conflicts = []
                # Calculate
                try:
                    while i < len(graph.index):
                        # Get the pathways of the node
                        node = graph.index[i]
                        node_pathways = {'activators': [], 'inhibitors': []}
                        [node_pathways['activators'].append(pathway) if pathway.activator
                         else node_pathways['inhibitors'].append(pathway) for pathway in pathways
                         if pathway.consequent == node]
                        # Solve the conflicts
                        pathways = list(filter(lambda x: x.consequent != node, pathways))
                        manager = ConflictsManager(activators=node_pathways['activators'],
                                                   inhibitors=node_pathways['inhibitors'],
                                                   priority_matrix=blosum,
                                                   network=initial_network,
                                                   conflicts=conflicts,
                                                   base_map=base_map,
                                                   graph=graph,
                                                   node=node)
                        pathways.extend(manager.get_solution())
                        pathways.sort(key=lambda x: x.consequent)
                        # INPUT validation:
                        # There cannot be any pathway with an input in the consequent
                        if len(pathways) != len(list(filter(lambda x: x.consequent not in inputs, pathways))):
                            raise ValueError
                        # Filtering of equivalent pathways
                        occurrences = []
                        codes = []
                        for p in range(0, len(pathways)):
                            path1 = pathways[p]
                            for q in range(0, len(pathways)):
                                path2 = pathways[q]
                                if path1.region_of_interest == path2.region_of_interest and\
                                        path1.consequent == path2.consequent:
                                    code = ''.join(path1.region_of_interest) + path1.consequent
                                    occurrences += [(p, code)]
                                    codes += [code] if code not in codes else []
                        occurrences = [list(filter(lambda x: x[1] == code, occurrences))[0] for code in codes]
                        pathways = [pathways[p[0]] for p in occurrences]
                        # Validation
                        if i == len(graph.index) - 1:
                            iter_count += 1
                            if iter_count >= max_iterations:
                                condition_result = False
                                break
                            if len(pathways) != n_pathways:
                                i = 0
                                n_pathways = len(pathways)
                                continue
                        i += 1
                except InputAlterationException:
                    # If the map of the INPUT is altered the simulation is not valid.
                    condition_result = False
                except ConvergenceException:
                    # If the resolution does not converge
                    condition_result = False
                # Check the condition and add the new result
                if condition_result:
                    # Apply all the pathways over the map
                    base_map.modificate_maps(pathways)
                    # Get the expression from the maps
                    expressions = base_map.get_expressions()
                    # Validation with the attractors
                    primes = FileExchange.bnet2primes('\n'.join(expressions))
                    stg = StateTransitionGraphs.primes2stg(primes, "synchronous")
                    steady, cyclic = Attractors.compute_attractors_tarjan(stg)
                    result = Result(network=initial_network,
                                    pathways=pathways,
                                    maps_set=base_map,
                                    conflicts=conflicts,
                                    simulation=m,
                                    iterations=iter_count,
                                    variant=None,
                                    expressions=expressions,
                                    attractors={'steady': steady, 'cyclic': cyclic})
                    condition = all([True if att in steady else False for att in attractors])
                    # Add the result
                    if condition:
                        results['accepted'].append(result)
                    else:
                        results['dismissed'].append(result)
            bar.update(m)
    print('Simulations completed')
    # Filtering for equivalent results
    codes = []
    final_results = []
    [(codes.append(result.code), final_results.append(result)) for result in results['accepted']
     if result.code not in codes]
    return final_results


def blosum_generator(graph):
    """
    DESCRIPTION:
    A function to randomly generate a priority matrix in the way of blosum matrices for alignments.
    :param graph: [pandas Dataframe] the graph to whom we aim to produce the matrix.
    :return: the randomly generated matrix.
    """
    labels = [''.join(sorted(it)) for sl
              in [combinations(graph.index, i) for i in range(1, len(graph.index) + 1)] for it in sl]
    content = [random.sample(range(1, 2 * len(labels) * len(labels) + 1), len(labels)) for row in labels]
    content = [[0 if i == j else content[i][j] for j in range(0, len(content[i]))] for i in range(0, len(content))]
    blosum = pd.DataFrame(data=content, index=labels, columns=labels)

    return blosum

def netFilter(networks):
    """
    DESCRIPTION:
    A function to let only the promising networks. The ones which can be useful.
    """
    # Auxiliary functions
    def simple_projection(nets, nodes):
        for net in nets:
            steadies = net['steady']
            network = net['network'].split('\n')
            positions = [i for i in range(0, len(network)) if network[i][0] in nodes]
            chunks = [it for sl in [[(pos, steady[pos]) for steady in steadies] for pos in positions] for it in sl]
            condition = all(
                [True if len(list(set(group[:]))) == 1 else False for group in
                 [[chunk[1] for chunk in chunks if chunk[0] == pos] for pos in positions]]
            )
            if condition:
                values = dict(set([(network[chunk[0]][0], chunk[1]) for chunk in chunks]))
                network = '\n'.join(
                    [network[i].split(',')[0] + ', ' + values[network[i].split(',')[0]]
                     if i in positions else network[i] for i in range(0, len(network))]
                )
                primes = FileExchange.bnet2primes(network)
                stg = StateTransitionGraphs.primes2stg(primes, "synchronous")
                steady, cyclic = Attractors.compute_attractors_tarjan(stg)
                net['steady'] = steady
                net['cyclic'] = cyclic
                yield net
    nodes = ['E', 'F']
    networks = list(simple_projection(nets=networks, nodes=nodes))
    print()


