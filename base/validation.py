#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import time
from functools import reduce
import itertools
import progressbar
from graphs.conflicts_theory import Pathway, ConflictsManager
from base.exceptions import InputAlterationException, ConvergenceException
from base.hibrids import Result
from base.ncbf import blosum_generator
from PyBoolNet import StateTransitionGraphs, Attractors, FileExchange
import multiprocessing as mp
import tqdm


class Validation:
    """
    DESCRIPTION:
    An object to standardize the different kinds of validation gathered in in ncbf.py.
    """

    def __init__(self, networks, nodes, inputs, kind='Tarjan', max_global_iterations=20, max_local_iterations=20,
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

    def validation(self, iteration):
        """
        DESCRIPTION:
        A method introduced to multiprocess the amount of validation to be performed.
        :param iteration: [tuple] with the number of simulations and the network object.
        :return: [Result] the result of the network creation and assessment. All the information.
        """
        m = iteration[0]
        network = iteration[1]
        # iterations_count = iteration_bar[0]
        # bar = iteration_bar[1][1]
        # Generate the priority matrix
        priority_matrix = blosum_generator(network.variant.data)
        condition_result = True  # set the condition to append the result of the simulation
        pathways = network.pathways[:]
        n_pathways = len(pathways)  # number of pathways
        # Control parameters
        i = 0
        iter_count = 0
        conflicts = []
        # Extra number of rounds to avoid unforeseen conflicts
        p = 0
        p_limit = 5
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
                # Redundant condition for unforeseen conflicts
                if i >= len(network.variant.data.index) and p <= p_limit:
                    i = 0
                    p += 1

        except InputAlterationException:
            # If the map of the INPUT is altered the simulation is not valid.
            condition_result = False
        except ConvergenceException:
            # If the resolution does not converge
            condition_result = False
        # Update progress
        # bar.update(iterations_count)
        # iterations_count += 1
        # Check the condition and add the new result
        if condition_result:
            # Apply all the pathways over the map
            network.original_map.modificate_maps(pathways)
            # Get the expression from the maps
            expressions = list(network.original_map.get_expressions())
            # Validation with the attractors
            primes = FileExchange.bnet2primes('\n'.join(expressions))
            stg = StateTransitionGraphs.primes2stg(primes, "synchronous")
            steady, cyclic = Attractors.compute_attractors_tarjan(stg)
            result = Result(network=network.structure,
                            pathways=pathways,
                            maps_set=network.original_map,
                            conflicts=conflicts,
                            simulation=m,
                            iterations=iter_count,
                            expressions=expressions,
                            attractors={'steady': steady, 'cyclic': cyclic},
                            variant=network.variant,
                            priority_matrix=priority_matrix,
                            roles_set=network.variant.roles)
            # Set if the result has been successful and return it
            result.accepted = all([True if att in steady else False for att in self.attractors])
            return result

    def third_validation(self, max_pathways=6400):
        """
        DESCRIPTION:
        The mechanism of validation to obtain the corrected maps with all the generated pathways according to the given,
        initial, map.
        :param method: [string] validation method for the attractors.
        :param max_pathways: [int] parameter to set the maximum number of pathways.
        """
        # Launch the simulations
        iters = list(itertools.product(range(0, self.simulations), self.networks, repeat=1))
        n_procs = mp.cpu_count()
        pool = mp.Pool(n_procs)
        jobs = pool.map_async(self.validation, iters)
        pool.close()
        results = list(filter(lambda x: x is not None, jobs.get()))
        # Delete repeated networks
        registry = []
        final_results = []
        for result in results:
            if result.expressions not in registry:
                registry.append(result.expressions)
                final_results.append(result)
        return final_results



    def group_results(self):
        """
        DESCRIPTION:
        This method is devised to group the networks according to its similarity to the one described through the
        attractors.
        :return: [dictionary] grouped results.
        """
        # Group non-cyclic attractors
        non_cyclic_results = [result for result in self.results if len(result.attractors['cyclic']) == 0]
        # With the same number of steady states and non cyclic
        same_number = [result for result in non_cyclic_results if len(result.attractors['steady']) == len(self.attractors)]
        # Group by number of steady attractors of the list
        at_least_one_of_all = []
        for result in non_cyclic_results:
            if any([True if att in result.attractors['steady'] else False for att in self.attractors]):
                at_least_one_of_all.append(result)
        sn_at_least_one_of_all = []
        for result in same_number:
            if any([True if att in result.attractors['steady'] else False for att in self.attractors]):
                sn_at_least_one_of_all.append(result)
        at_least_two_of_all = []
        for result in at_least_one_of_all:
            if len([att for att in self.attractors if att in result.attractors['steady']]) > 1:
                at_least_two_of_all.append(result)
        sn_at_least_two_of_all = []
        for result in sn_at_least_one_of_all:
            if len([att for att in self.attractors if att in result.attractors['steady']]) > 1:
                sn_at_least_two_of_all.append(result)
        at_least_three_of_all = []
        for result in at_least_two_of_all:
            if len([att for att in self.attractors if att in result.attractors['steady']]) > 2:
                at_least_three_of_all.append(result)
        sn_at_least_three_of_all = []
        for result in sn_at_least_two_of_all:
            if len([att for att in self.attractors if att in result.attractors['steady']]) > 2:
                sn_at_least_three_of_all.append(result)
        at_least_four_of_all = []
        for result in at_least_three_of_all:
            if len([att for att in self.attractors if att in result.attractors['steady']]) > 3:
                at_least_four_of_all.append(result)
        sn_at_least_four_of_all = []
        for result in sn_at_least_three_of_all:
            if len([att for att in self.attractors if att in result.attractors['steady']]) > 3:
                sn_at_least_four_of_all.append(result)
        # Organize the response
        grouped_results = {
            'total_results': self.results,
            'non_cyclic_results': non_cyclic_results,
            'same_number': same_number,
            'one_of_all': at_least_one_of_all,
            'same_number_with_at_least_one_of_all': sn_at_least_one_of_all,
            'two_of_all': at_least_two_of_all,
            'same_number_with_at_least_two_of_all': sn_at_least_two_of_all,
            'three_of_all': at_least_three_of_all,
            'same_number_and_three_of_all': sn_at_least_three_of_all,
            'four_of_all': at_least_four_of_all,
            'same_number_and_four_of_all': sn_at_least_four_of_all
        }
        return grouped_results

    def execute_validation(self):
        """
        DESCRIPTION:
        A method to perform the validation according to the different validation methods (Tarjan o STP).
        :return: [list] the results of the validation.
        """
        self.results = list(self.third_validation())
        # Group results by similarity
        return self.group_results()

