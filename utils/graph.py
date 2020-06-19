import copy
from functools import reduce
import pandas as pd
import uuid
import itertools
from alive_progress import alive_bar
from utils.conflicts_theory import KMap, Pathway
from utils.ncbf import ncbfCalc, networksCalc
from utils.validation import Validation


class Graph:
    """
    DESCRIPTION:
    An object to represent an instance of the represented graph. With all its possible manifestations.
    """

    # Methods
    def __init__(self, initial_data, attractors=None, filter_kernel=None, imposed_roles_sets=None, simulations=20,
                 variants_limit=None, max_global_iterations=20, max_local_iterations=20):
        """
        DESCRIPTION:
        Builder of the object graph.
        :param initial_data: [pandas DataFrame] data representing the initial data upon which the graph is built.
        :param attractors: [list] strings containing the attractors for the validation.
        :param filter_kernel: [dictionary] all the data to perform the filtration.
        :param imposed_roles_sets: [list] roles sets imposed externally, not directly inferred from the graph.
        :param simulations: [int] number of simulations to be performed in every network.
        :param variants_limit: [int] limit to the number of variants to be generated.
        :param max_global_iterations: [int] parameter to set the maximum number of iterations in every network.
        :param max_local_iterations: [int] parameter to set the maximum number of iterations for the conflicts solver.
        """
        # Set the attributes
        self.initial_data = initial_data
        self.networks_codes = []
        self.filter = Filter(**filter_kernel) if filter_kernel is not None else Filter()
        self.nodes = self.get_nodes()
        self.space = self.get_space()
        self.attractors = attractors
        # Impose the given roles sets
        self.imposed_roles_sets = []
        if imposed_roles_sets is not None:
            [list(map(lambda x: x.append(str(x[-2]) + str(x[-1])), roles_set)) for roles_set in imposed_roles_sets]
            self.imposed_roles_sets = imposed_roles_sets if imposed_roles_sets is not None else []
        # Establish the relationships between nodes
        self.set_adjustment()
        self.inputs = self.get_inputs()
        # Obtain the roles over the basis structure
        self.roles_combinations = self.get_roles_combinations()
        # Set the variants of the network
        print('Generating variants')
        self.variants = self.get_variants(limit=variants_limit)
        print('Variants generation completed')
        print('Searching for networks')
        self.networks = self.get_networks()
        print('Networks search completed')
        print('Launch the validation of the networks')
        self.validation = Validation(networks=self.networks, nodes=[node.name for node in self.get_nodes()],
                                     inputs=self.inputs, attractors=self.attractors, simulations=simulations,
                                     max_global_iterations=max_global_iterations,
                                     max_local_iterations=max_local_iterations)
        # Obtain the results of the validation
        self.results = self.validation.get_results()
        print()

    def __str__(self):
        """
        DESCRIPTION:
        String method of the object
        :return: [string] a readable representation of the object.
        """
        return f'Graph nodes {self.nodes}'

    def get_nodes(self):
        """
        DESCRIPTION:
        A method to generate at once all the nodes of the map.
        :return: [list] nodes of the graph.
        """
        if 'nodes' in dir(self):
            # If the nodes exist, return them. They are to be adjusted.
            nodes = self.nodes
        else:
            # If not, generate them. They are not adjusted.
            nodes = [Node(name=node) for node in self.initial_data.index]
        return nodes

    def get_node(self, name=None, id=None):
        """
        DESCRIPTION:
        A method to return node by ID or name.
        :param id: [string] id of the searched node.
        :param name: [string] name of the searched node.
        :return: [Node] searched node.
        """
        node = list(filter(lambda x: x.name == name or x.id == id, self.nodes))[0]
        return node

    def get_roles_combinations(self):
        """
        DESCRIPTION:
        A method to get all the possible combinations with the associated roles in the nodes.
        :return: [list] sets of combinations.
        """
        # Auxiliary functions
        def aux_fun1(group):
            sources = [row[0] for row in group]
            value = True if len(set(sources)) == len(group) else False
            return value

        def aux_fun2(group):
            sources = [row[1] for row in group]
            value = True if len(set(sources)) == len(group) else False
            return value

        def aux_fun3(tables):
            codes = [[]] * len(tables)
            for i in range(0, len(tables)):
                if len(tables[i]) == 1:
                    yield [tables[i]]
                    continue
                set_rows = list(aux_fun4(codes, tables[i], i))
                yield list(filter(lambda x: x != 0, set_rows))

        def aux_fun4(codes, row, index):
            for j in range(0, len(row)):
                if j == 0:
                    codes[index] = [0] * len(row)
                new_row = sorted(row[j], key=lambda x: x[1])
                if str(new_row) not in codes[index]:
                    codes[index][j] = str(new_row)
                    yield new_row

        # Assess if we have already defined the attribute or not
        if 'roles_combinations' not in dir(self):
            # Make all the virtual networks with the roles tables
            tables = [node.roles_table.groupby('node') for node in self.get_nodes()]
            tables = [[it for sl in [table.get_group(group).values.tolist() for group in table.groups] for it in sl]
                      for table in tables]
            tables = [
                list(
                    filter(
                        aux_fun2,
                        itertools.product(tables[i],
                                          repeat=len(self.get_nodes()[i].activators + self.get_nodes()[i].inhibitors)))
                )
                if len(tables[i]) > 1 else tables[i] for i in range(0, len(tables))
            ]
            tables = list(aux_fun3(tables))
            combinations = tables[0]
            for i in range(1, len(tables)):
                combinations = [list(map(lambda x: x + [str(x[2]) + str(x[3])], [it for sl in group for it in sl]))
                                for group in itertools.product(combinations, tables[i], repeat=1)]
            # At the end we impose the extra roles sets
            self.roles_combinations = list(map(lambda x: [it[0:5] for it in x], combinations)) + self.imposed_roles_sets

        return self.roles_combinations

    def get_variants(self, limit=None):
        """
        DESCRIPTION:
        A method to obtain all the variants of the graph.
        :param limit: [int] number of variants to which the generation can be limited.
        :return: [list] associated variants.
        """
        if 'variants' not in dir(self):
            variants = list(self.variants_generator(limit))
            self.variants = variants
        return self.variants

    def get_inputs(self):
        """
        DESCRIPTION
        A method to get all the inputs of the graph.
        """
        if 'inputs' not in dir(self):
            self.inputs = [node.name for node in list(filter(lambda x: x.get_input(), self.get_nodes()))]
        return self.inputs

    def get_networks(self):
        """
        DESCRIPTION:
        A method to get all the different object networks associated with the graph through the variants.
        :return: [list] associated networks.
        """
        if 'networks' not in dir(self):
            # Get total number of networks
            self.networks = list(self.networks_generator())
        return self.networks

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

        nodes = [node.name for node in self.get_nodes()]
        if 'space' not in dir(self):
            combinations = [bin(i).split('b')[1] if len(bin(i).split('b')[1]) == len(nodes) else
                            ''.join(str_gen(len(nodes) - len(bin(i).split('b')[1]))) + bin(i).split('b')[1]
                            for i in range(0, 2 ** len(nodes))]
            self.space = [{nodes[i]: int(c[i]) for i in range(0, len(c))} for c in combinations]

        return self.space

    def set_adjustment(self):
        """
        DESCRIPTION:
        A method to establish the relationship between nodes.
        """
        # Auxiliary functions
        def aux_fun(node_name):
            node = self.get_node(name=node_name)
            [node.set_associated_node(ass_node, True) for ass_node in self.get_nodes()
             if ass_node.name in self.initial_data.loc[node_name, 'activators']]
            [node.set_associated_node(ass_node, False) for ass_node in self.get_nodes()
             if ass_node.name in self.initial_data.loc[node_name, 'inhibitors']]
        # Set the relationships between nodes
        list(map(aux_fun, self.initial_data.index))
        # Set the roles table of every node
        [node.set_roles_table() for node in self.get_nodes()]

    def variants_generator(self, limit=None):
        """
        DESCRIPTION:
        A method to generate all the variants associated with the graph.
        :param limit: [int] number of variants to which the generation can be limited.
        """
        inputs_names = self.get_inputs()
        roles_combinations = self.get_roles_combinations()
        roles_combinations = roles_combinations[0:limit] if limit is not None else roles_combinations
        with alive_bar(len(roles_combinations)) as bar:
            for i in range(0, len(roles_combinations)):
                if self.filter.clean(roles_set=roles_combinations[i]):
                    yield Variant(roles=roles_combinations[i], initial_data=self.initial_data, inputs=inputs_names,
                                  space=self.get_space())
                bar()

    def networks_generator(self):
        """
        DESCRIPTION:
        A method to generate all the networks associated with the variants in the graph.
        """
        with alive_bar() as bar:
            for variant in self.get_variants():
                for network in variant.get_networks():
                    if self.filter.clean(structure=network):
                        # The copy is needed to avoid the recursive application over all pathways
                        yield Network(structure=network, variant=copy.deepcopy(variant))
                    bar()


class Node:
    """
    DESCRIPTION:
    A class to represent elements of the network to be represented.
    """

    # Methods
    def __init__(self, name):
        """
        DESCRIPTION:
        Builder of the class.
        :param name: [string] name of the node.
        """
        self.id = f'nd{uuid.uuid1()}'
        self.name = name
        self.activators = []
        self.inhibitors = []
        self.roles_table = None

    def __str__(self):
        """
        DESCRIPTION:
        String method of the object
        :return: [string] a readable representation of the object.
        """
        return f'Node: ID {self.id} Name {self.name}'

    def get_input(self):
        """
        DESCRIPTION:
        :return: [boolean] indicator if the node is an input or not.
        """
        return True if not self.activators and not self.inhibitors else False

    def set_associated_node(self, node, activator):
        """
        DESCRIPTION:
        A method to set both the activator and inhibitor associated nodes of the given node.
        :param node: [Node] node to be introduced.
        :param activator: [boolean] nature of the node with respect this.
        """
        if activator:
            self.activators.append(node)
        else:
            self.inhibitors.append(node)

    def set_roles_table(self):
        """
        DESCRIPTION:
        A method to build the table which sets the variants in behavior of every node associated with this.
        Notes:
            - node: the node with the behaviour.
            - canalizing: the intended value of the variable.
            - canalized: the intended value of the function.
        """
        if not self.get_input():
            data = []
            [data.extend([[self.name, node.name, 1, 1], [self.name, node.name, 0, 0]]) for node in self.activators]
            [data.extend([[self.name, node.name, 0, 1], [self.name, node.name, 1, 0]]) for node in self.inhibitors]
        else:
            data = [[self.name, self.name, 1, 1]]
        self.roles_table = pd.DataFrame(data=data, columns=['source', 'node', 'canalizing', 'canalized'])


class Variant:
    """
    DESCRIPTION:
    An object to set a modification of the graph depending the roles of activators and inhibitors.
    """

    # Methods
    def __init__(self, roles, initial_data, inputs, space):
        """
        DESCRIPTION:
        Builder of the class.
        :param roles: [list] representation of the role given to every node.
        :param initial_data: [pandas DataFrame] initial representation of the graph.
        :param inputs: [list] names of the nodes that act as inputs.
        :param space: [list] dictionaries representing the space in which the expressions are to be assessed
        """
        self.roles = roles
        self.space = space
        self.initial_data = initial_data
        self.inputs = inputs
        self.data = self.get_data()
        self.pathways = list(self.get_initial_pathways())
        self.paths = self.get_paths()
        self.networks = self.get_networks()
        self.roles_code = self.get_roles_code()

    def __str__(self):
        """
        DESCRIPTION:
        String method of the object
        :return: [string] a readable representation of the object.
        """
        return f'Variant roles code: {self.get_roles_code()}'

    def get_data(self):
        """
        DESCRIPTION:
        A method to modify the initial data given by the graph accordingly with the roles.
        :return: [pandas DataFrame] modified representation of the graph.
        """
        # Parameters
        tags = {'act': ['11', '01'], 'inh': ['00', '10']}
        data = pd.DataFrame(index=self.initial_data.index, columns=self.initial_data.columns)
        # Generate and return the modified data
        if 'data' not in dir(self):
            for node in data.index:
                if node in self.inputs:
                    data.at[node, :] = [[''], ['']]
                    continue
                act = list(map(lambda x: x[1], filter(lambda x: x[0] == node and x[-1] in tags['act'], self.roles)))
                inh = list(map(lambda x: x[1], filter(lambda x: x[0] == node and x[-1] in tags['inh'], self.roles)))
                data.loc[node, 'activators'] = act if act else ['']
                data.loc[node, 'inhibitors'] = inh if inh else ['']
        else:
            data = self.data
        return data

    def get_paths(self):
        """
        DESCRIPTION:
        A method to obtain the basic paths of the network. The first step to launch the tree algorithm to obtain the
        networks.
        :return: [pandas Series] a path o every kind.
        """
        if 'paths' not in dir(self):
            tags = ['activators', 'inhibitors']
            paths = ncbfCalc(data=self.data, tags=tags)
        else:
            paths = self.paths
        return paths

    def get_networks(self):
        """
        DESCRIPTION:
        A method to get all the networks given the set of paths in the graph.
        :return: [list] all networks.
        """
        if 'networks' not in dir(self):
            networks = list(networksCalc(self.paths))
        else:
            networks = self.networks
        return networks

    def get_initial_pathways(self):
        """
        DESCRIPTION:
        A method to obtain the initial set of pathways given a variant. Beware! It is assumed that if a node is not in
        the activators, it is to be in the inhibitors and vice versa.
        """
        for i in range(0, len(self.data.index)):
            for el in list(reduce(lambda x, y: x + y, self.data.values.tolist()[i])):
                if el != '':
                    if el in self.data['activators'][self.data.index[i]]:
                        yield Pathway(antecedent=el, consequent=self.data.index[i], activator=True, space=self.space)
                    else:
                        yield Pathway(antecedent=el, consequent=self.data.index[i], activator=False, space=self.space)

    def get_roles_code(self):
        """
        DESCRIPTION:
        A method to, given a roles set, generate its code.
        """
        if 'roles_code' not in dir(self):
            self.roles.sort(key=lambda x: x[0] + x[1])
            self.roles_code = ''.join(list(map(lambda x: x[0] + x[1] + str(x[2]) + str(x[3]), self.roles)))
        return self.roles_code


class Network:
    """
    DESCRIPTION:
    An object to represent every network to be validated.
    """

    # Methods
    def __init__(self, structure, variant):
        """
        DESCRIPTION:
        Builder method of the object.
        :param structure: [list] structure of the network to be taken by the algorithm.
        :param variant: [Variant] object with the variant of the graph
        """
        self.structure = structure
        self.variant = variant
        self.map = self.get_map()
        self.code = self.get_code()
        self.pathways = self.get_pathways()
        self.set_canalizing_values()

    def __str__(self):
        """
        DESCRIPTION:
        Method to obtain a readable representation of the object.
        :return: [string] a readable representation of the object.
        """
        exprs = '$'.join([path.expression for path in self.pathways])
        return '|'.join([' '.join(step) for step in self.structure if step != 'INPUT']) + ' @ ' + exprs

    def get_map(self):
        """
        DESCRIPTION:
        A method to build the object map associated with every network.
        :return: [KMap] table of truth of the network.
        """
        if 'map' not in dir(self):
            self.map = KMap(network=self.structure, graph=self.variant.data, roles_set=self.variant.roles,
                            inputs=self.variant.inputs, space=self.variant.space)
        return self.map

    def get_code(self):
        """
        DESCRIPTION:
        A method to obtain or generate the code which characterizes the kind to which the network belongs.
        """
        def code_maker(item):
            return f'${"".join(item[1])}:{item[0]}$'

        if 'code' not in dir(self):
            self.code = ''.join(map(code_maker, self.map.get_support(hex_flag=True).items()))
        return self.code

    def get_pathways(self):
        """
        DESCRIPTION:
        A method to generate the network pathways imposing the given role over the variant pathways.
        """
        if 'pathways' not in dir(self):
            # Initial pathways
            self.pathways = self.variant.pathways
        return self.pathways

    def set_canalizing_values(self):
        """
        DESCRIPTION:
        A method to impose over the expressions of the pathways the canalizing values.
        """
        groups = filter(lambda x: x[1][0] == x[0].consequent and x[1][1] == x[0].antecedent,
                        itertools.product(self.get_pathways(), self.variant.roles, repeat=1))
        self.pathways = list(map(impose_roles, groups))


class Filter:
    """
    DESCRIPTION:
    An object to control the filtration over the networks of the graph through the control in the generation of variants
    and their associated networks.
    """

    # Methods
    def __init__(self, roles_sets=None, structures=None):
        """
        DESCRIPTION:
        Builder of the class.
        :param roles_sets: [list] with the roles sets allowed.
        :param structures: [list] with the structures allowed.
        """
        self.roles_sets_codes = list(map(self.roles_code_maker, roles_sets)) if roles_sets is not None else []
        self.structures_codes = list(map(self.structures_code_maker, structures)) if structures is not None else []

    def roles_code_maker(self, roles_set):
        """
        DESCRIPTION:
        A method to, given a roles set, generate its code.
        :param roles_set: [list] roles of the different nodes.
        """
        roles_set.sort(key=lambda x: x[0] + x[1])
        value = ''.join(list(map(lambda x: x[0] + x[1] + str(x[2]) + str(x[3]), roles_set)))
        return value

    def structures_code_maker(self, structure):
        """
        DESCRIPTION:
        A method to, given a roles set, generate its code.
        :param structure: [list] structure of the network.
        """
        return '@'.join(
            [
                '$'.join([''.join(sorted(step)) if len(step) > 1 else step for step in path])
                if path != 'INPUT'
                else path for path in structure
            ]
        )

    def clean(self, **kwargs):
        """
        DESCRIPTION:
        A method to perform the validation given by the filter.
        """
        condition = True
        # Filtering by roles set
        if 'roles_set' in kwargs.keys() and self.roles_sets_codes:
            kwargs['roles_set'].sort(key=lambda x: x[0])
            condition = condition and self.roles_code_maker(kwargs['roles_set']) in self.roles_sets_codes
        # Filtering by structure
        if 'structure' in kwargs.keys() and self.structures_codes:
            code = self.structures_code_maker(kwargs['structure'])
            condition = condition and code in self.structures_codes
        return condition


def impose_roles(group):
    """
    DESCRIPTION:
    A function to efficiently modify the expressions of the set of pathways.
    :param group: [tuple] the pathway and its associated role.
    """
    group[0].set_role(group[1])
    return group[0]
