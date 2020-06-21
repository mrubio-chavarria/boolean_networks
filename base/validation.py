#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from functools import reduce
from pytictoc import TicToc
from alive_progress import alive_bar
from graphs.conflicts_theory import Pathway, ConflictsManager
from base.exceptions import InputAlterationException, ConvergenceException
from base.hibrids import Result
from base.ncbf import blosum_generator
from PyBoolNet import StateTransitionGraphs, Attractors, FileExchange


class Validation:
    """
    DESCRIPTION:
    An object to standardize the different kinds of validation gathered in in ncbf.py.
    """

    def __init__(self, networks, nodes, inputs, kind='III', max_global_iterations=20, max_local_iterations=20,
                 attractors=None, simulations=20):
        """
        DESCRIPTION:
        Constructor of the object Validation.
        :param kind: [string] code indicating the kind of validation to be performed.
        :param networks: [list] networks to be validated.
        :param nodes: [list] names of the nodes of the network.
        :param inputs: [list] nodes that are inputs.
        :param attractors: [list] strings containing the attractors for the validation.
        :param simulations: [int] number of simulations to be performed in the conflicts resolution.
        :param max_global_iterations: [int] parameter to set the maximum number of iterations all around the network.
        :param max_local_iterations: [int] parameter to set the maximum number of iterations in conflicts solver.
        """
        self.type = kind
        self.networks = networks
        self.nodes = nodes
        self.inputs = inputs
        self.simulations = simulations
        self.max_global_iterations = max_global_iterations
        self.max_local_iterations = max_local_iterations
        self.attractors = attractors
        self.space = self.get_space()
        # Obtain the results of the validation
        self.results = self.execute_validation()

    def get_space(self):
        """
        DESCRIPTION:
        A method to establish the boolean space in which the variants are going to be evaluated.
        :return: [list] dictionaries representing the space.
        """
        # Auxiliary functions
        def str_gen(n):
            for i in range(0, n):
                yield '0'

        # Generate the space
        if 'space' not in dir(self):
            combinations = [bin(i).split('b')[1] if len(bin(i).split('b')[1]) == len(self.nodes) else
                            ''.join(str_gen(len(self.nodes) - len(bin(i).split('b')[1]))) + bin(i).split('b')[1]
                            for i in range(0, 2 ** len(self.nodes))]
            return [{self.nodes[i]: int(c[i]) for i in range(0, len(c))} for c in combinations]
        else:
            return self.space

    def get_initial_pathways(self, variant):
        """
        DESCRIPTION:
        A method to obtain the initial set of pathways given a variant. Beware! It is assumed that if a node is not in
        the activators, it is to be in the inhibitors and vice versa.
        :param variant: [Variant] object variant with all the information.
        """
        for i in range(0, len(self.nodes)):
            for el in list(reduce(lambda x, y: x + y, variant.data.values.tolist()[i])):
                if el != '':
                    if el in variant.data['activators'][self.nodes[i]]:
                        yield Pathway(antecedent=el, consequent=self.nodes[i], activator=True, space=self.get_space())
                    else:
                        yield Pathway(antecedent=el, consequent=self.nodes[i], activator=False, space=self.get_space())

    def third_validation(self, max_pathways=6400):
        """
        DESCRIPTION:
        The mechanism of validation to obtain the corrected maps with all the generated pathways according to the given,
        initial, map.
        :param max_pathways: [int] parameter to set the maximum number of pathways.
        """
        # Launch the simulations
        with alive_bar(self.simulations*len(self.networks)) as bar:
            for m in range(0, self.simulations):
                for network in self.networks:
                    # Generate the priority matrix
                    priority_matrix = blosum_generator(network.variant.data)
                    condition_result = True  # set the condition to append the result of the simulation
                    pathways = network.pathways[:]
                    n_pathways = len(pathways)  # number of pathways
                    # Control parameters
                    i = 0
                    iter_count = 0
                    conflicts = []
                    # Calculate
                    try:
                        while i < len(network.variant.data.index):
                            # Get the pathways of the node
                            node = network.variant.data.index[i]
                            node_pathways = {'activators': [], 'inhibitors': []}
                            [node_pathways['activators'].append(pathway) if pathway.activator
                             else node_pathways['inhibitors'].append(pathway) for pathway in pathways
                             if pathway.consequent == node]
                            # Solve the conflicts
                            pathways = list(filter(lambda x: x.consequent != node, pathways))
                            manager = ConflictsManager(activators=node_pathways['activators'],
                                                       inhibitors=node_pathways['inhibitors'],
                                                       priority_matrix=priority_matrix,
                                                       network=network.structure,
                                                       conflicts=conflicts,
                                                       base_map=network.map,
                                                       graph=network.variant.data,
                                                       algorithm='I',
                                                       max_iterations=self.max_local_iterations,
                                                       node=node)
                            pathways.extend(manager.get_solution())
                            pathways.sort(key=lambda x: x.consequent)
                            # INPUT validation:
                            # There cannot be any pathway with an input in the consequent
                            if len(pathways) != len(list(filter(lambda x: x.consequent not in self.inputs, pathways))):
                                raise InputAlterationException
                            # Filtering of equivalent pathways
                            occurrences = []
                            codes = []
                            for p in range(0, len(pathways)):
                                path1 = pathways[p]
                                for q in range(0, len(pathways)):
                                    path2 = pathways[q]
                                    if path1.region_of_interest == path2.region_of_interest and \
                                            path1.consequent == path2.consequent:
                                        code = ''.join(path1.region_of_interest) + path1.consequent
                                        occurrences += [(p, code)]
                                        codes += [code] if code not in codes else []
                            occurrences = [list(filter(lambda x: x[1] == code, occurrences))[0] for code in codes]
                            pathways = [pathways[p[0]] for p in occurrences]
                            # Validation
                            if i == len(network.variant.data.index) - 1:
                                iter_count += 1
                                if iter_count >= self.max_global_iterations:
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
                    # Update progress
                    bar()
                    # Check the condition and add the new result
                    if condition_result:
                        # Apply all the pathways over the map
                        network.map.modificate_maps(pathways)
                        # Get the expression from the maps
                        t = TicToc()
                        expressions = list(network.map.get_expressions())
                        # Validation with the attractors
                        primes = FileExchange.bnet2primes('\n'.join(expressions))
                        stg = StateTransitionGraphs.primes2stg(primes, "synchronous")
                        steady, cyclic = Attractors.compute_attractors_tarjan(stg)
                        result = Result(network=network.structure,
                                        pathways=pathways,
                                        maps_set=network.map,
                                        conflicts=conflicts,
                                        simulation=m,
                                        iterations=iter_count,
                                        expressions=expressions,
                                        attractors={'steady': steady, 'cyclic': cyclic},
                                        variant=network.variant)
                        # Set if the result has been successful and return it
                        result.accepted = all([True if att in steady else False for att in self.attractors])
                        yield result

    def group_results(self):
        """
        DESCRIPTION:
        This method is devised to group the networks according to its similarity to the one described through the
        attractors.
        :return: [dictionary] grouped results.
        """
        # Group by cyclic and steady attractors presence
        non_cyclic_results = [result for result in self.results if len(result.attractors['cyclic']) == 0]
        # With the same number of steady states and non cyclic
        steadies = [result for result in non_cyclic_results if len(result.attractors['steady']) == len(self.attractors)]
        # Group by number of steady attractors of the list
        one_of_all = [result for result in steadies
                      if len([att for att in self.attractors if att in result.attractors]) == 1]
        two_of_all = [result for result in steadies
                      if len([att for att in self.attractors if att in result.attractors]) == 2]
        three_of_all = [result for result in steadies if result.accepted]
        # Organize the response
        grouped_results = {
            'total_results': self.results,
            'non_cyclic_results': non_cyclic_results,
            'same_number': steadies,
            'one_of_all': one_of_all,
            'two_of_all': two_of_all,
            'three_of_all': three_of_all
        }
        return grouped_results

    def execute_validation(self):
        """
        DESCRIPTION:
        A method to perform the validation according to the different validation methods.
        :return: [list] the results of the validation.
        """

        # Select between the different methods of validation
        if self.type == 'I':
            # First method of validation
            self.results = []
        elif self.type == 'II':
            # Second method of validation
            self.results = []
        if self.type == 'III':
            # Third method of validation
            self.results = list(self.third_validation())
        # Group results by similarity
        return self.group_results()

