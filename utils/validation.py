
import uuid
import pytictoc
from functools import reduce
from utils.conflicts_theory import Pathway, KMap
from utils.ncbf import blosum_generator


class Validation:
    """
    DESCRIPTION:
    An object to standardize the different kinds of validation gathered in in ncbf.py.
    """

    def __init__(self, variants, nodes, inputs, kind='III'):
        """
        DESCRIPTION:
        Builder of the object Validation.
        :param kind: [string] code indicating the kind of validation to be performed.
        :param variants: [list] variants to be validated.
        :param nodes: [list] names of the nodes of the network.
        :param inputs: [list] nodes that are inputs.
        """
        self.id = f'vd{uuid.uuid1()}'
        self.type = kind
        self.variants = variants
        self.nodes = nodes
        self.inputs = inputs
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
            yield [Pathway(antecedent=el, consequent=self.nodes[i], activator=True, space=self.get_space())
                   if el in variant.data['activators'][self.nodes[i]] else
                   Pathway(antecedent=el, consequent=self.nodes[i], activator=False, space=self.get_space())
                   for el in list(reduce(lambda x, y: x + y, variant.data.values.tolist()[i])) if el != '']

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
                    yield 3

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


