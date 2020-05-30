import re
import sys
import itertools
import pandas as pd
from utils.utils import Term
from utils.Kmap import Minterms


class Pathway:
    """
    DESCRIPTION:
    An object to represent the pathway with which we are working.
    """

    # Methods
    def __init__(self, antecedent, consequent, activator, space=None, expression=None):
        """
        DESCRIPTION:
        Constructor of the object.
        :param antecedent: [string] left side of the equation which describes the pathway. Condition.
        :param consequent: [string] right side of the equation which describes the pathway. Result of the condition.
        :param activator: [boolean] Variable to store if the relation among sides if activatory or inhibitory.
        :param space: [list] list of dicts with all variables combinations which form space in which the function is
        defined.
        """
        self.antecedent = ''.join(sorted(antecedent))
        self.consequent = consequent
        self.activator = activator
        self.expression = self.set_antecedent_expression_from_graph(antecedent) if expression is None else expression
        if space is None:
            self.map = None
            self.region_of_interest = None
        else:
            self.set_map(space)

    def __str__(self):
        return self.antecedent + ' --> ' + self.consequent + '\n Activator: ' + str(self.activator)

    def set_expression(self, psi, graph, variables):
        """
        DESCRIPTION:
        On the contrary to one below this method is designed to correct the expression. If it is required it will make
        another pathway in order to keep the expression coherent.
        :param psi: [list] minterms produced by karnaugh maps simplification.
        :param graph: [pandas DataFrame] the graph with which we are working.
        :param variables: [list] dicts representing the space, variables combinations, in which the expression is to be
        assessed.
        :return: [list] pathways resulting from the modification of the expression.
        """
        # Auxiliary functions
        def minterm_checker(factors):
            conditions = []
            for factor_1 in factors[:]:
                if factor_1[0] == '!':
                    continue
                conditions.extend([True if (factor_1 not in factor_2) or (factor_1 in factor_2 and '!' not in factor_2)
                                   else False for factor_2 in factors[:]])
            return all(conditions)

        # Correct the expression and make all the combinations in the left side of the pathway
        psi = [[graph.index[i] if word[i] == '0' else '!' + graph.index[i]
                for i in range(0, len(word)) if word[i] != '*'] for word in psi]
        # psi = [['!B', '!C'], ['!A', '!B'], ['!D', 'C']]
        if len(psi) != 1:
            psi = [[[[(f1, f2) for f2 in psi[j]] for f1 in psi[i]] for j in range(i+1, len(psi))] for i in range(0, len(psi)-1)]
            psi = [list(it1) for sl1 in [it2 for sl2 in [it3 for sl3 in psi for it3 in sl3] for it2 in sl2] for it1 in sl1]
        new_psi = []
        [new_psi.append(it) for sl in [[[var, self.expression] for var in group] for group in psi] for it in sl
         if minterm_checker(it) if it not in new_psi]
        psi = [list(set(term)) for term in new_psi]

        # Make all the necessary pathways
        r = re.compile('\w')
        new_pathways = []
        for group in psi:
            expression = '&'.join(group)
            antecedent = ''.join(r.findall(expression))
            new_pathways.append(Pathway(antecedent=antecedent, consequent=self.consequent, activator=self.activator,
                                        space=variables, expression=expression))

        self.antecedent = ''.join(r.findall(self.expression))
        # Update the region of interest
        self.set_map(variables)
        return new_pathways

    def set_antecedent_expression_from_graph(self, original):
        """
        DESCRIPTION:
        A method to develop a expression of the antecedent from the graph. In other words, this method is not designed
        to introduced the corrections of the zones.
        :param original: [list] al the nodes involved in the condition.
        :return: [string] the logical expression in boolnet format.
        """
        expression = '&'.join(list(original))
        return expression

    def eval_expression(self, variables):
        """
        DESCRIPTION:
        A method to assess whether the expression of the pathway meets the condition or not.
        :param variables: [dict] variables and their values.
        :return: [boolean] the value of the expression.
        """
        # Auxiliary functions
        def local_eval(minterm, variables):
            condition = lambda x: variables[x[0]] == 1 if x[0] != '!' else not variables[x[1]] == 1
            if '!' == minterm[0]:
                minterm = minterm[1::]
                response = not all(map(condition, minterm))
            else:
                response = all(map(condition, minterm))
            return response

        values = []
        r = re.compile(r'\!\([^\)]+\)|\w(?![^(]*\))|[!](?=\w)')
        minterms = [item for item in r.findall(self.expression)]
        for i in range(0, len(minterms)):
            if minterms[i] == '!':
                minterms[i+1] = '!' + minterms[i+1]

        minterms = list(filter(lambda x: x != '!', minterms))
        values = [local_eval(minterm, variables) for minterm in minterms]
        value = all(values)
        return value

    def set_map(self, variables_set):
        """
        DESCRIPTION:
        A method that, given a set of variables calculates its Karnaugh map in an abstract manner.
        :param variables_set: [list] list of dicts with the variables combinations to be tested.
        :return: [dict] combinations and their result to the pathway.
        :return: [list] strings drawing the to which the pathway has been designed.
        """
        self.map = {''.join([str(var) for var in vs_set.values()]): self.eval_expression(vs_set)
                    for vs_set in variables_set}
        self.region_of_interest = [key for key in self.map.keys() if self.map[key]]


class Conflict:
    """
    DESCRIPTION:
    An object to represent the conflict between two pathways.
    """

    # Methods
    def __init__(self, first_pathway, second_pathway, priority_matrix, psi, graph):
        """
        DESCRIPTION:
        Constructor of the object.
        :param first_pathway: [pathway] The first pathway interfering in the conflict.
        :param second_pathway: [pathway] The second pathway interfering in the conflict.
        :param priority_matrix: [pandas DataFrame] Matrix with the scores to solve the conflict.
        :param psi: [set] region of the space in which the first and second pathways overlap.
        :param graph: [pandas DataFrame] description of the whole set of pathways from which these two comes.
        """
        self.first_pathway = first_pathway
        self.second_pathway = second_pathway
        try:
            self.priority = self.set_priority(priority_matrix) # Pathway with the lowest priority.
        except:
            print(first_pathway)
            print(second_pathway)
        self.psi = psi
        self.graph = graph

    def __str__(self):
        return f'Conflict:\n {self.first_pathway}\n vs\n {self.second_pathway}'

    def set_priority(self, priority_matrix):
        """
        DESCRIPTION:
        Priority of A vs B corresponds with A (row) against B (column). In this field we put the pathway with the lowest
        priority.
        :param priority_matrix: [pandas DataFrame] Matrix with the scores to determine which pathway has the lowest
        priority and, consequently, it is to be modified.
        """
        firsts_priority = priority_matrix.loc[self.first_pathway.antecedent, self.second_pathway.antecedent]
        seconds_priority = priority_matrix.loc[self.second_pathway.antecedent, self.first_pathway.antecedent]
        if firsts_priority < seconds_priority:
            self.priority = self.first_pathway
        else:
            self.priority = self.second_pathway
        return self.priority

    def solve(self, initial_network):
        """
        DESCRIPTION:
        A method to solve the conflict through the modification of the pathway with lowest priority. It too creates
        new pathways.
        :param initial_network: [list] the structure in which the graph nodes are distributed.
        """
        # Auxiliary functions
        def str_gen(n):
            for i in range(0, n):
                yield '0'
        combinations = [bin(i).split('b')[1] if len(bin(i).split('b')[1]) == len(self.graph.index) else
                        ''.join(str_gen(len(self.graph.index) - len(bin(i).split('b')[1]))) + bin(i).split('b')[1]
                        for i in range(0, 2 ** len(self.graph.index))]
        variables = [{self.graph.index[i]: int(comb[i]) for i in range(0, len(comb))} for comb in combinations]
        # Simplify psi
        psi = [Term(value) for value in list(self.psi)]
        psi = Minterms(psi)
        psi.simplify()
        psi = [str(term) for term in psi.result]
        original_psi = [f"{'&'.join([self.graph.index[i] if word[i] == '1' else '!' + self.graph.index[i] for i in range(0, len(word)) if word[i] != '*'])}" for word in psi]
        original_psi = '&'.join(original_psi)
        # Correct the pathway and generate pathways resulting from the correction
        new_pathways = self.priority.set_expression(psi, self.graph, variables)
        new_pathways = new_pathways if new_pathways is not None else []
        high_pathway = [path for path in [self.first_pathway, self.second_pathway] if path != self.priority][0]
        # Generate the basis map
        kmap = KMap(initial_network, self.graph)
        # Impose the solution over the map
        high_pathway.expression = high_pathway.antecedent
        minterms = [''.join([str(value) for value in variable.values()]) for variable in variables
                    if high_pathway.eval_expression(variable)]
        kmap.maps.at[high_pathway.consequent, minterms] = high_pathway.activator
        positions = kmap.maps.loc[high_pathway.consequent, :]
        positions = positions[positions != high_pathway.activator].index
        minterms = [Term(minterm) for minterm in positions]
        minterms = Minterms(minterms)
        minterms.simplify()
        s = [str(term) for term in minterms.result]
        s_nb = [f"{'&'.join([self.graph.index[i] if word[i] == '1' else '!' + self.graph.index[i] for i in range(0, len(word)) if word[i] != '*'])}"
                for word in s]
        # Generate new pathways
        r = re.compile('\w')
        for minterm in s_nb:
            for var in minterm.split('&'):
                if var[0] == '!':
                    var = var[1]
                    activator = False
                else:
                    var = var[0]
                    activator = True
                antecedent = r.findall(original_psi)  # The antecedent is supposed to be built upon the letter only
                new_pathways.append(Pathway(antecedent=antecedent, consequent=var, activator=activator,
                                            space=variables))

        return new_pathways


class KMap:
    """
    DESCRIPTION:
    An object to represent the map of a given network.
    """

    # Methods
    def __init__(self, network, graph):
        """
        DESCRIPTION:
        Constructor of the object.
        :param network: [list] the structure in which the graph nodes are distributed.
        :param graph: [pandas DataFrame] description of the whole set of pathways from which the network comes.
        """
        self.network = network
        self.graph = graph
        self.maps = self.set_maps()

    def eval_expression(self, variables_set, expressions):
        """
        DESCRIPTION:
        An eval method to assess the expression of a set of expressions related with the nodes. Therefore, these
        expressions show a specific structure. It is not a general-purpose eval method.
        :param variables_set: [list] dicts with the combinations of variables to be assessed.
        :param expressions: [list] boolean functions to be assessed with the structure given by Murrugarra2013 in Boolnet
        notation.
        :return: [pandas DataFrame] all the functions associated to the nodes with the result for every variables
        combination.
        """
        # Auxiliary functions
        def local_eval(minterm, variables):
            condition = lambda x: variables[x[0]] == 1 if x[0] != '!' else not variables[x[1]] == 1
            if '!' == minterm[0]:
                minterm = minterm[2:-1].split('&')
                value = not all(map(condition, minterm))
            else:
                minterm = minterm.split('&')
                value = all(map(condition, minterm))
            return value

        final_values = []
        for expression in expressions:
            for variables in variables_set:
                # Manage the parentheses
                if '(' in expression:
                    r = re.compile('(?<=\()(.*)(?=\))')
                    parentheses_chunks = r.findall(expression)
                    parentheses_chunks = self.eval_expression([variables], parentheses_chunks)
                    r = re.compile('(?:^|\))([^()]*)(?:\(|$)')
                    non_parentheses_chunks = [i for i in r.findall(expression) if i != '']
                    # expr = [0] * (len(parentheses_chunks) + len(non_parentheses_chunks))
                    # Introduce correction for abandoned parentheses
                    np_content = [chunk if chunk[-1] != '!' else (chunk[0:-2], chunk[-1]) for chunk in
                                  non_parentheses_chunks]
                    new_chunks = []
                    [new_chunks.append(chunk) if type(chunk) == int else new_chunks.extend(chunk) for chunk in np_content]
                    new_chunks = [chunk for chunk in new_chunks if chunk != '']
                    expr = new_chunks + parentheses_chunks
                    for i in range(0, len(expr)):
                        if type(expr[i]) == str and expr[i] != '!':
                            factors = expr[i].split('&')
                            expr[i] = all([variables[factor[0]] == 1 if factor[0] != '!' else
                                           not variables[factor[1]] == 1 for factor in factors])
                    for i in range(0, len(expr)):
                        if expr[i] == '!':
                            expr[i + 1] = not expr[i + 1]
                    expr = [expr[i] for i in range(0, len(expr)) if expr[i] != '!']
                    final_values.append(all(expr))
                else:
                    factors = expression.split('&')
                    final_values.append(all([variables[factor[0]] == 1 if factor[0] != '!' else
                                             not variables[factor[1]] == 1 for factor in factors]))
        return final_values

    def set_maps(self):
        """
        DESCRIPTION:
        A method to calculate the response of the function in the whole space of the boolean expressions of the network.
        :return: [pandas DataFrame] the maps for all the nodes.
        """
        # Auxiliary functions
        def str_gen(n):
            for i in range(0, n):
                yield '0'

        combinations = [bin(i).split('b')[1] if len(bin(i).split('b')[1]) == len(self.graph.index) else
                        ''.join(str_gen(len(self.graph.index) - len(bin(i).split('b')[1]))) + bin(i).split('b')[1]
                        for i in range(0, 2 ** len(self.graph.index))]
        variables = [{self.graph.index[i]: int(comb[i]) for i in range(0, len(comb))} for comb in combinations]
        # Build the expression
        network = []
        for i in range(0, len(self.network)):
            # Set the structure of an INPUT node
            if self.network[i] == 'INPUT':
                self.network[i] = f'{self.graph.index[i]}'
                continue
            # Iterate through the layers of a node
            load = ''
            steps = range(len(self.network[i]) - 1, -1, -1)
            m_step = min(steps)
            first_letter = list(self.network[i][m_step])[0]
            for j in steps:
                layer = list(reversed(self.network[i][j]))
                # Calculate the whole load
                if len(layer) == 1:
                    level = '!' + layer[0]
                else:
                    level = layer[0]
                    for letter in layer[1::]:
                        level = '!' + letter + '&!' + level
                # Add inner load
                level = level + '&' + load if load != '' else level
                # Add the header
                if j != m_step or j == m_step and first_letter in self.graph.loc[self.graph.index[i], 'activators']:
                    load = '!(' + level + ')'
                else:
                    load = level
            network.append(load)
        # Build the map
        values = self.eval_expression(variables, network)
        values = [values[int(i*len(values)/len(self.graph.index)):int((i+1)*len(values)/len(self.graph.index))]
                  for i in range(0, len(self.graph.index))]
        index = self.graph.index
        columns = combinations
        # Change and send the map
        kmaps = pd.DataFrame(data=values, index=index, columns=columns)
        self.maps = kmaps
        return kmaps



