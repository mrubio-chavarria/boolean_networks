
import pandas as pd
import uuid
import itertools


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
        Constructor of the object graph.
        :param initial_data: [pandas DataFrame] data representing the initial data upon which the graph is built.
        """
        # Set the attributes
        self.id = f'gr{uuid.uuid1()}'
        self.initial_data = initial_data
        self.nodes = self.get_nodes()
        # Establish relationships between nodes
        self.set_adjustment()
        # Build related networks
        self.networks = self.get_networks()
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
        :return: nodes of the graph.
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

    def get_networks(self):
        """
        DESCRIPTION:
        A method to get all the networks with the associated roles in its nodes.
        :return: [list] networks with all the attributes.
        """
        # Auxiliary functions
        def aux_funI(group):
            sources = [row[0] for row in group]
            value = True if len(set(sources)) == len(group) else False
            return value

        def aux_funII(group):
            sources = [row[1] for row in group]
            value = True if len(set(sources)) == len(group) else False
            return value

        def aux_funIII(tables):
            codes = [[]] * len(tables)
            for i in range(0, len(tables)):
                if len(tables[i]) == 1:
                    yield [tables[i]]
                    continue
                set_rows = list(aux_funIV(codes, tables[i], i))
                yield list(filter(lambda x: x != 0, set_rows))

        def aux_funIV(codes, row, index):
            for j in range(0, len(row)):
                if j == 0:
                    codes[index] = [0] * len(row)
                new_row = sorted(row[j], key=lambda x: x[1])
                if str(new_row) not in codes[index]:
                    codes[index][j] = str(new_row)
                    yield new_row

        tables = [node.roles_table.groupby('node') for node in self.get_nodes()]
        tables = [[it for sl in [table.get_group(group).values.tolist() for group in table.groups] for it in sl]
                  for table in tables]
        tables = [
            list(
                filter(
                    aux_funII,
                    itertools.product(tables[i],
                                      repeat=len(self.get_nodes()[i].activators + self.get_nodes()[i].inhibitors)))
            )
            if len(tables[i]) > 1 else tables[i] for i in range(0, len(tables))
        ]
        tables = list(aux_funIII(tables))
        combinations = tables[0]
        for i in range(1, len(tables)):
            combinations = list(itertools.product(combinations, tables[i], repeat=2))
            print()

        return 3

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
        Constructor of the class.
        :param name: [string] name of the node.
        :param activators: [list] activator nodes of this one.
        :param inhibitors: [list] inhibitor nodes of this one.
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


