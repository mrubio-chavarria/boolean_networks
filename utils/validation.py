
import uuid
import pytictoc
from functools import reduce
from utils.conflicts_theory import Pathway, KMap, ConflictsManager
from utils.hibrids import Result
from utils.ncbf import blosum_generator
from PyBoolNet import StateTransitionGraphs, Attractors, FileExchange


class Validation:
    """
    DESCRIPTION:
    An object to standardize the different kinds of validation gathered in in ncbf.py.
    """

    def __init__(self, variants, nodes, inputs, kind='III', attractors=None):
        """
        DESCRIPTION:
        Builder of the object Validation.
        :param kind: [string] code indicating the kind of validation to be performed.
        :param variants: [list] variants to be validated.
        :param nodes: [list] names of the nodes of the network.
        :param inputs: [list] nodes that are inputs.
        :param attractors: [list] strings containing the attractors for the validation.
        """
        self.id = f'vd{uuid.uuid1()}'
        self.type = kind
        self.variants = variants
        self.nodes = nodes
        self.inputs = inputs
        self.attractors = attractors
        self.space = self.get_space()
        self.results = self.execute_validation()

    def __str__(self):
        """
        DESCRIPTION:
        String method of the object
        :return: [string] a readable representation of the object.
        """
        return f'Validation: ID {self.id}'

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

    def third_validation(self, simulations=20, max_iterations=2000, max_pathways=6400):
        """
        DESCRIPTION:
        The mechanism of validation to obtain the corrected maps with all the generated pathways according to the given,
        initial, map.
        :param simulations: [int] parameter to set the number of simulations to perform.
        :param max_iterations: [int] parameter to set the maximum number of iterations all around the network.
        :param max_pathways: [int] parameter to set the maximum number of pathways.
        """
        # Validate each variant
        for variant in self.variants:
            # Obtain the initial set of pathways from the graph
            initial_pathways = list(self.get_initial_pathways(variant=variant))
            # Launch the simulations
            for m in range(0, simulations):
                # Generate the matrix with the priorities for the simulation
                blosum = blosum_generator(variant.data)
                net_counter = -1
                for network in variant.get_networks():
                    net_counter += 1
                    condition_result = True  # set the condition to append the result of the simulation
                    pathways = initial_pathways[:]
                    n_pathways = len(pathways)  # number of pathways
                    # Generate the basis map
                    base_map = KMap(network=network, graph=variant.data, roles_set=variant.roles, inputs=self.inputs)
                    # Control parameters
                    i = 0
                    iter_count = 0
                    conflicts = []
                    # Calculate
                    try:
                        while i < len(variant.data.index):
                            # Get the pathways of the node
                            node = variant.data.index[i]
                            node_pathways = {'activators': [], 'inhibitors': []}
                            [node_pathways['activators'].append(pathway) if pathway.activator
                             else node_pathways['inhibitors'].append(pathway) for pathway in pathways
                             if pathway.consequent == node]
                            # Solve the conflicts
                            pathways = list(filter(lambda x: x.consequent != node, pathways))
                            manager = ConflictsManager(activators=node_pathways['activators'],
                                                       inhibitors=node_pathways['inhibitors'],
                                                       priority_matrix=blosum,
                                                       network=network,
                                                       conflicts=conflicts,
                                                       base_map=base_map,
                                                       graph=variant.data,
                                                       node=node)
                            pathways.extend(manager.get_solution())
                            pathways.sort(key=lambda x: x.consequent)
                            # INPUT validation:
                            # There cannot be any pathway with an input in the consequent
                            if len(pathways) != len(list(filter(lambda x: x.consequent not in self.inputs, pathways))):
                                raise ValueError
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
                            if i == len(variant.data.index) - 1:
                                iter_count += 1
                                if iter_count >= max_iterations:
                                    condition_result = False
                                    break
                                if len(pathways) != n_pathways:
                                    i = 0
                                    n_pathways = len(pathways)
                                    continue
                            i += 1
                    except ValueError:
                        # If the map of the INPUT is altered the simulation is not valid. Likewise, if we enter in a
                        # cyclic resolution of the conflicts, the simulation is not valid.
                        condition_result = False
                    except IndexError:
                        # The same situation of the INPUT, at the time of selecting destiny nodes.
                        condition_result = False
                    except RecursionError:
                        # Another situation in which the error algorithm takes too many iterations.
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
                        result = Result(network=network,
                                        pathways=pathways,
                                        maps_set=base_map,
                                        conflicts=conflicts,
                                        simulation=m,
                                        iterations=iter_count,
                                        expressions=expressions,
                                        attractors={'steady': steady, 'cyclic': cyclic})
                        # Set if the result has been successful and return it
                        result.accepted = all([True if att in steady else False for att in self.attractors])
                        yield result

    def execute_validation(self):
        """
        DESCRIPTION:
        A method to perform the validation according to the different validation methods.
        :return: [list] the results of the validation.
        """
        # Select between the different methods of validation
        if self.type == 'I':
            # First method of validation
            pass
        elif self.type == 'II':
            # Second method of validation
            pass
        if self.type == 'III':
            # Third method of validation
            results = list(self.third_validation())
        return results


