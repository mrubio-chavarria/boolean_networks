
import pandas as pd
import uuid
import itertools
from utils.ncbf import ncbfCalc, networksCalc
import progressbar


class Graph:
    """
    DESCRIPTION:
    An object to represent an instance of the represented graph. With all its possible manifestations.
    """
    # Attributes

    # Methods
    def __init__(self, initial_data):
        """
        DESCRIPTION:
        Builder of the object graph.
        :param initial_data: [pandas DataFrame] data representing the initial data upon which the graph is built.
        """
        # Set the attributes
        self.id = f'gr{uuid.uuid1()}'
        self.initial_data = initial_data
        self.nodes = self.get_nodes()
        # Establish the relationships between nodes
        self.set_adjustment()
        self.inputs = self.get_inputs()
        # Obtain the roles over the basis structure
        self.roles_combinations = self.get_roles_combinations()
        # Set the variants of the network
        print('Start variants generation')
        self.variants = self.get_variants()
        print('Variants generation completed')
        print()

    def __str__(self):
        """
        DESCRIPTION:
        String method of the object
        :return: [string] a readable representation of the object.
        """
        return f'Graph: ID {self.id} Nodes {self.nodes}'

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
            combinations = list(map(lambda x: [it[0:5] for it in x], combinations))
        else:
            combinations = self.roles_combinations

        return combinations

    def get_variants(self):
        """
        DESCRIPTION:
        A method to obtain all the variants of the graph.
        :return: [list] associated variants.
        """
        if 'variants' not in dir(self):
            variants = list(self.variants_generator())
            self.variants = variants
        return self.variants

    def get_inputs(self):
        """
        DESCRIPTION
        A method to get all the inputs of the graph.
        """
        if 'inputs' not in dir(self):
            self.inputs = list(filter(lambda x: x.get_input(), self.get_nodes()))
        return self.inputs

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

    def variants_generator(self):
        """
        DESCRIPTION:
        A method to generate all the variants associated with the graph.
        """
        inputs_names = [node.name for node in self.get_inputs()]
        roles_combinations = self.get_roles_combinations()
        with progressbar.ProgressBar(max_value=len(roles_combinations)) as bar:
            for i in range(0, len(roles_combinations)):
                yield Variant(roles=roles_combinations[i], initial_data=self.initial_data, inputs=inputs_names)
                bar.update(i)


class Node:
    """
    DESCRIPTION:
    A class to represent elements of the network to be represented.
    """
    # Attributes

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
    # Attributes

    # Methods
    def __init__(self, roles, initial_data, inputs):
        """
        DESCRIPTION:
        Builder of the class.
        :param roles: [list] representation of the role given to every node.
        :param initial_data: [pandas DataFrame] initial representation of the graph.
        :param inputs: [list] names of the nodes that act as inputs.
        """
        self.id = f'vr{uuid.uuid1()}'
        self.roles = roles
        self.initial_data = initial_data
        self.inputs = inputs
        self.data = self.get_data()
        self.paths = self.get_paths()
        self.networks = self.get_networks()

    def __str__(self):
        """
        DESCRIPTION:
        String method of the object
        :return: [string] a readable representation of the object.
        """
        return f'Variant: ID {self.id}'

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
