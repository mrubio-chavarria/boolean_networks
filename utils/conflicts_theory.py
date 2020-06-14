import re
import sys
import itertools
import pandas as pd
import time
from utils.exceptions import InputAlterationException
from utils.utils import Term
from utils.Kmap import Minterms
import random
import sympy
from sympy import Symbol
from sympy.logic import SOPform
from uuid import uuid4


class ConflictsManager:
    """
    DESCRIPTION:
    An object to handle the conflicts resolution in a row.
    """

    # Attributes

    # Methods
    def __init__(self, activators, inhibitors, priority_matrix, graph, conflicts, network, base_map, node):
        """
        DESCRIPTION:
        Builder of the class.
        :param activators: [list] pathways to be considered as activators.
        :param inhibitors: [list] pathways to be considered as inhibitors.
        :param priority_matrix: [pandas DataFrame] matrix with the priority to solve the conflicts.
        :param graph: [pandas DataFrame] graph with the nodes of the network.
        :param network: [list] structure of the network.
        :param base_map: [pandas DataFrame] representation of the truth table of all nodes.
        :param node: [string] node to which belong the conflicts solved.
        :param conflicts: [list] group of conflicts to be extended.
        """
        self.activators = activators
        self.inhibitors = inhibitors
        self.priority_matrix = priority_matrix
        self.graph = graph
        self.network = network
        self.base_map = base_map
        self.node = node
        self.extra_pathways = []
        self.registry = {}
        self.conflicts = conflicts
        self.following_group = {'activators': [], 'inhibitors': []}
        self.algorithm_counter = 0
        if self.activators != [] and self.inhibitors != []:
            self.set_registry()
            # list(self.conflicts_resolution_algorithm())
            self.launch_algorithm()

    def set_registry(self, duple=None):
        """
        DESCRIPTION:
        Method to establish the registry at any time given the accumulated data.
        :param duple: [tuple] pair to be registered.
        """
        # Initial case
        if duple is None:
            if not self.registry:
                self.registry = {key: [] for key in [path.id for path in self.activators + self.inhibitors]}
            else:
                [self.registry.update({key: []}) for key in [path.id for path in self.activators + self.inhibitors]
                 if key not in self.registry.keys()]
        elif duple is not None:
            if duple[0].id in self.registry.keys():
                self.registry[duple[0].id].append(duple[1])
            else:
                self.registry[duple[0].id] = duple[1]
            if duple[1].id in self.registry.keys():
                self.registry[duple[1].id].append(duple[0])
            else:
                self.registry[duple[1].id] = duple[0]

    def launch_algorithm(self):
        """
        DESCRIPTION:
        The method which holds the recurrent algorithm to solve the conflicts, provided a set of pathways.
        """
        # Initial parameters
        working_group = []
        current_inhibitors = []
        pathways_limit = 100
        counter_limit = 1000
        # Delete repeated pathways in activators, inhibitors and extra_pathways
        codes = []
        activators = []
        inhibitors = []
        extra = []
        [(codes.append(pathway.code), activators.append(pathway))
         if pathway.activator else (codes.append(pathway.code), inhibitors.append(pathway))
         for pathway in self.activators + self.inhibitors if pathway.code not in codes]
        [(codes.append(pathway.code), extra.append(pathway)) for pathway in self.extra_pathways
         if pathway.code not in codes]
        self.activators = activators
        self.inhibitors = inhibitors
        self.extra_pathways = extra
        # Add the following group of the last iteration
        [self.activators.append(pathway) if pathway.activator else self.inhibitors.append(pathway)
         for pathway in self.following_group['inhibitors'] + self.following_group['activators']]
        self.following_group['activators'] = []
        self.following_group['inhibitors'] = []
        # Set the working group
        for act in self.activators:
            inhs = [inh for inh in self.inhibitors if inh not in self.registry[act.id] and inh not in current_inhibitors]
            inh = random.choice(inhs) if inhs else None
            if inh is None:
                self.following_group['activators'].append(act)
                continue
            current_inhibitors.append(inh)
            working_group.append((act, inh))
        if not self.activators:  # it wont happen
            self.following_group['inhibitors'].extend(self.inhibitors)
        # Set the non-paired inhibitors in the next round
        [self.following_group['inhibitors'].append(inh) for inh in self.inhibitors if inh not in current_inhibitors]
        # Solve the conflicts
        new_pathways = []
        change_pathways = [False]*len(working_group)
        for i in range(0, len(working_group)):
            # Set evidence in registry
            self.set_registry(working_group[i])
            # Calculate conflicts
            psi = set(working_group[i][0].region_of_interest) & set(working_group[i][1].region_of_interest)
            if len(psi) > 0:
                conflict = Conflict(working_group[i][0], working_group[i][1], priority_matrix=self.priority_matrix,
                                    psi=psi, graph=self.graph)
                self.conflicts.append(conflict)
                new_pathways.extend(conflict.solve2(self.network, self.base_map))
                change_pathways[i] = True
            else:
                new_pathways.extend([working_group[i][0], working_group[i][1]])
        # Redistribute the new pathways
        if any(change_pathways):
            self.activators = list(filter(lambda x: x.consequent == self.node and x.activator, new_pathways))
            self.inhibitors = list(filter(lambda x: x.consequent == self.node and not x.activator, new_pathways))
            # Add the extra pathways
            self.add_extra_pathways(new_pathways)
            # Register the new pathways
            self.set_registry()
            # New pathways, launch again the algorithm
            self.algorithm_counter += 1
            if self.algorithm_counter > counter_limit:
                raise RecursionError
            self.launch_algorithm()
        else:
            completed_activators = [False]*len(self.activators)
            if len(self.activators + self.inhibitors) > pathways_limit:
                # If we enter into a cyclic generation of pathways we are to stop. TO BE IMPROVED.
                raise ValueError
            for i in range(0, len(self.activators)):
                unregistered_inhibitors = list(
                    filter(lambda x: x not in self.registry[self.activators[i].id], self.inhibitors)
                )
                if not unregistered_inhibitors:
                    completed_activators[i] = True
            if not all(completed_activators):
                # There pathways to compare, launch the algorithm again
                self.algorithm_counter += 1
                if self.algorithm_counter > counter_limit:
                    raise RecursionError
                self.launch_algorithm()

    def conflicts_resolution_algorithm(self, activators=None, inhibitors=None, following=None):
        """
        DESCRIPTION:
        The improved version of the algorithm above.
        :return: [list] final set of pathways.
        """
        if activators is not None and inhibitors is not None:
            # Case in which we are inside the nested iteration
            pairs = pd.DataFrame(columns=['activators', 'inhibitors'])
            for activator in activators:
                pass
        else:
            # Case in which the method is called for the first time
            # Get the initial parameters
            activators = self.activators[:]
            inhibitors = self.inhibitors[:]
            # Solve
            yield from self.conflicts_resolution_algorithm(activators=activators, inhibitors=inhibitors)


    def get_solution(self):
        """
        DESCRIPTION:
        A method to return the response, the final set of pathways: extra, activators and inhibitors.
        :return: [list] the final pathways.
        """
        # Compound the set
        pathways = self.extra_pathways + self.activators + self.inhibitors
        # Filtering for the repeated and equivalent pathways
        # TO BE IMPLEMENTED IF NECESSARY
        return pathways

    def add_extra_pathways(self, local_pathways):
        """
        DESCRIPTION:
        A method to filter the addition of new pathways in every addition. And to avoid the repeated ones.
        :params pathways: [list] pathways to be added in extra pathways.
        """
        self.extra_pathways.extend(local_pathways)
        self.extra_pathways = list(filter(lambda x: x.consequent != self.node, self.extra_pathways))
        # Filtering for repeated extra pathways
        codes = []
        extra_pathways = []
        [(codes.append(pathway.code), extra_pathways.append(pathway)) for pathway in self.extra_pathways
         if pathway.code not in codes]


class Pathway:
    """
    DESCRIPTION:
    An object to represent the pathway with which we are working.
    """

    # Methods
    def __init__(self, antecedent, consequent, activator, space=None, expression=None, canalizing_value=None,
                 canalized_value=None):
        """
        DESCRIPTION:
        Builder of the object.
        :param antecedent: [string] left side of the equation which describes the pathway. Condition.
        :param consequent: [string] right side of the equation which describes the pathway. Result of the condition.
        :param activator: [boolean] Variable to store if the relation among sides if activatory or inhibitory.
        :param space: [list] list of dicts with all variables combinations which form space in which the function is
        defined.
        :param canalizing_value: [int] value to be manifested by the antecedent.
        :param canalized_value: [int] value to be provoked in the consequent.
        """
        self.id = str(uuid4())
        self.code = ''
        self.antecedent = ''.join(sorted(antecedent))
        self.canalizing_value = canalizing_value
        self.canalized_value = canalized_value
        self.consequent = consequent
        self.activator = activator
        self.expression = self.set_antecedent_expression_from_graph(antecedent) if expression is None else expression
        if space is not None:
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
            condition = True
            new_factors = []
            [new_factors.append(it) for sl in [factor.split('&') for factor in factors]
             for it in sl if it not in new_factors]
            for factor_1 in new_factors[:]:
                for factor_2 in new_factors[:]:
                    if factor_1 == '!' + factor_2 or factor_2 == '!' + factor_1:
                        condition = False
            return condition

        # Correct the expression and make all the combinations in the left side of the pathway
        psi = [[graph.index[i] if word[i] == '0' else '!' + graph.index[i]
                for i in range(0, len(word)) if word[i] != '*'] for word in psi]
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

    def set_map(self, variables_set):
        """
        DESCRIPTION:
        A method that, given a set of variables calculates its Karnaugh map in an abstract manner.
        :param variables_set: [list] list of dicts with the variables combinations to be tested.
        :return: [dict] combinations and their result to the pathway.
        :return: [list] strings drawing the to which the pathway has been designed.
        :return: [string] code to indicate the type to which the pathway belongs.
        """
        self.map = {''.join([str(var) for var in vs_set.values()]): self.eval_expression(vs_set)
                    for vs_set in variables_set}
        self.region_of_interest = list(sorted([key for key in self.map.keys() if self.map[key]]))
        self.code = f'${"".join(self.region_of_interest)}:{self.consequent}$'

    def set_role(self, role):
        """
        DESCRIPTION:
        A method to impose upon the expression the canalizing value. The canalizing value has been taken into account
        before.
        """
        factors = self.expression.split('&')
        for i in range(0, len(factors)):
            if role[1] in factors[i]:
                if len(factors[i]) == 2 and role[2] == 0:
                    factors[i] = factors[i][-1]
                elif len(factors[i]) == 1 and role[2] == 1:
                    factors[i] = '!' + factors[i]
                    factors[i] = '!' + factors[i]
        self.expression = '&'.join(factors)

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

        r = re.compile(r'\!\([^\)]+\)|\w(?![^(]*\))|[!](?=\w)')
        minterms = [item for item in r.findall(self.expression)]
        for i in range(0, len(minterms)):
            if minterms[i] == '!':
                minterms[i+1] = '!' + minterms[i+1]

        minterms = list(filter(lambda x: x != '!', minterms))
        values = [local_eval(minterm, variables) for minterm in minterms]
        value = all(values)
        return value


class Conflict:
    """
    DESCRIPTION:
    An object to represent the conflict between two pathways.
    """

    # Methods
    def __init__(self, first_pathway, second_pathway, priority_matrix, psi, graph):
        """
        DESCRIPTION:
        Builder of the object.
        :param first_pathway: [pathway] The first pathway interfering in the conflict.
        :param second_pathway: [pathway] The second pathway interfering in the conflict.
        :param priority_matrix: [pandas DataFrame] Matrix with the scores to solve the conflict.
        :param psi: [set] region of the space in which the first and second pathways overlap.
        :param graph: [pandas DataFrame] description of the whole set of pathways from which these two comes.
        """
        self.first_pathway = first_pathway
        self.second_pathway = second_pathway
        self.priority = self.set_priority(priority_matrix) # Pathway with the lowest priority.
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
        :return: [Pathway] the pathway object with the lowest priority.
        """
        firsts_priority = priority_matrix.loc[self.first_pathway.antecedent, self.second_pathway.antecedent]
        seconds_priority = priority_matrix.loc[self.second_pathway.antecedent, self.first_pathway.antecedent]
        if firsts_priority < seconds_priority:
            self.priority = self.first_pathway
        else:
            self.priority = self.second_pathway
        return self.priority

    def solve(self, initial_network, base_map):
        """
        DESCRIPTION:
        A method to solve the conflict through the modification of the pathway with lowest priority. It too creates
        new pathways.
        :param initial_network: [list] the structure in which the graph nodes are distributed.
        :param base_map: [KMap] the object with the information of the map upon which all the conflicts are being
        applied.
        :return: [dict] the whole set of pathways to be deleted and to be added to the rest.
        """
        # Auxiliary functions
        def str_gen(n):
            for i in range(0, n):
                yield '0'
        combinations = [bin(i).split('b')[1] if len(bin(i).split('b')[1]) == len(self.graph.index) else
                        ''.join(str_gen(len(self.graph.index) - len(bin(i).split('b')[1]))) + bin(i).split('b')[1]
                        for i in range(0, 2 ** len(self.graph.index))]
        variables = [{self.graph.index[i]: int(comb[i]) for i in range(0, len(comb))} for comb in combinations]
        # Simplify psi: THIS SECTION IS TO BE UPDATED TO THE NEW SIMPLIFICATION METHOD
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
        # Impose the solution over the map
        # new_pathways.extend(high_pathway.set_expression(psi, self.graph, variables))
        # high_pathway.expression = high_pathway.antecedent
        minterms = [''.join([str(value) for value in variable.values()]) for variable in variables
                    if high_pathway.eval_expression(variable)]
        base_map.maps.at[high_pathway.consequent, minterms] = high_pathway.activator
        positions = base_map.maps.loc[high_pathway.consequent, :]
        positions = positions[positions != high_pathway.activator].index
        minterms = [Term(minterm) for minterm in positions]
        minterms = Minterms(minterms)
        minterms.simplify()
        s = [str(term) for term in minterms.result]
        # Select the most advantageous term
        # TO BE IMPLEMENTED. By this time it is only selected a random term.
        word = random.choice(s)
        minterm = f"{'&'.join([self.graph.index[i] if word[i] == '1' else '!' + self.graph.index[i] for i in range(0, len(word)) if word[i] != '*'])}"
        # Generate new pathways
        r = re.compile('\w')
        for var in minterm.split('&'):
            if var[0] == '!':
                var = var[1]
                activator = False
            else:
                var = var[0]
                activator = True
            antecedent = r.findall(original_psi)  # The antecedent is supposed to be built upon the letter only
            new_pathways.append(Pathway(antecedent=antecedent, consequent=var, activator=activator, space=variables))
        # Create the response
        response = {
            'delete': self.priority,
            'add': new_pathways
        }
        return response

    def solve2(self, initial_network, base_map):
        """
        DESCRIPTION:
        A method to solve the conflict through the modification of the pathway with lowest priority. It too creates
        new pathways.
        :param initial_network: [list] the structure in which the graph nodes are distributed.
        :param base_map: [KMap] the object with the information of the map upon which all the conflicts are being
        applied.
        :return: [dict] the whole set of pathways to be deleted and to be added to the rest.
        """
        # Auxiliary functions
        def str_gen(n):
            for i in range(0, n):
                yield '0'
        combinations = [bin(i).split('b')[1] if len(bin(i).split('b')[1]) == len(self.graph.index) else
                        ''.join(str_gen(len(self.graph.index) - len(bin(i).split('b')[1]))) + bin(i).split('b')[1]
                        for i in range(0, 2 ** len(self.graph.index))]
        variables = [{self.graph.index[i]: int(comb[i]) for i in range(0, len(comb))} for comb in combinations]
        # Obtain the inputs
        inputs = [self.graph.index[i] for i in range(0, len(initial_network)) if initial_network[i] == 'INPUT']
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
        # Impose the solution over the map
        # new_pathways.extend(high_pathway.set_expression(psi, self.graph, variables))
        # high_pathway.expression = high_pathway.antecedent
        minterms = [''.join([str(value) for value in variable.values()]) for variable in variables
                    if high_pathway.eval_expression(variable)]
        base_map.maps.at[high_pathway.consequent, minterms] = high_pathway.activator
        positions = base_map.maps.loc[high_pathway.consequent, :]
        minterms = [int(position[0], 2) for position in positions[positions != high_pathway.activator].index]
        simplified_minterms = str(simplify(list(self.graph.index), minterms)).replace(' ', '').replace('~', '!').\
            split('|')
        # Select the most advantageous term
        # TO BE IMPLEMENTED. By this time it is only selected a random term.
        r = re.compile('\w')
        valid_minterms = [minterm for minterm in simplified_minterms
                          if not any([True if letter in inputs else False for letter in r.findall(minterm)])]
        if not valid_minterms:
            # If there are no valid terms, the algorithm has terminated
            raise InputAlterationException
        minterm = random.choice(valid_minterms)
        # Generate new pathways
        for var in minterm.split('&'):
            if var[0] == '!':
                var = var[1]
                activator = False
            else:
                var = var[0]
                activator = True
            antecedent = r.findall(original_psi)  # The antecedent is supposed to be built upon the letter only
            new_pathways.append(Pathway(antecedent=antecedent, consequent=var, activator=activator, space=variables))
        # Send the pathways
        new_pathways.append(high_pathway)
        # Modificate the map with the new pathways
        base_map.modificate_maps(new_pathways)
        return new_pathways


class KMap:
    """
    DESCRIPTION:
    An object to represent the map of a given network.
    """

    # Methods
    def __init__(self, network, graph, roles_set=None, inputs=None, space=None):
        """
        DESCRIPTION:
        Builder of the object.
        :param network: [list] the structure in which the graph nodes are distributed.
        :param graph: [pandas DataFrame] description of the whole set of pathways from which the network comes.
        :param roles_set: [list] the set of roles (canalizing and canalized value) for every node in the expressions.
        :param inputs: [list] names of the nodes of the network which are inputs.
        :param space: [list] dictionaries representing the space in which the expressions are to be assessed.
        """
        self.network = network
        self.graph = graph
        self.inputs = inputs
        self.roles_set = roles_set
        self.space = self.get_space() if space is None else space
        self.maps = self.set_maps()

    def __str__(self):
        """
        DESCRIPTION:
        A method to return a readable representation of the object.
        :return: [string] readable representation of the object.
        """
        maps = str(self.maps)
        roles = str(self.roles_set)
        return f'Maps: {maps} Roles: {roles}'

    def eval_expression(self, variables_set, expressions):
        """
        DESCRIPTION:
        An eval method to assess every expression related to the nodes. Therefore, these
        expressions show a specific structure. It is not a general-purpose eval method.
        :param variables_set: [generator] dicts with the combinations of variables to be assessed.
        :param expressions: [list] boolean functions to be assessed with the structure given by Murrugarra2013 in
        Boolnet notation.
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

    def modificate_maps(self, pathways):
        """
        DESCRIPTION:
        A method to modificate the set of map through the direct introduction of a set of pathways
        :param pathways: [list] pathways to be implemented over the map. Without any conflict. Directly.
        """
        def aux_fun(pathway):
            self.maps.loc[pathway.consequent, list(pathway.region_of_interest)] = pathway.activator
        list(map(aux_fun, pathways))

    def get_expressions(self):
        """
        DESCRIPTION:
        A function to extract the expressions stored in the maps set.
        :return: the expressions in BoolNet format of the nodes encoded in the set of maps.
        """
        # Obtain the letters
        nodes = self.maps.index
        # Simplification with KMaps
        for node in nodes:
            minterms = [it[0] for it in filter(lambda x: x[1], enumerate(self.maps.loc[node, :].tolist()))]
            sop = simplify(list(self.graph.index), minterms)
            if sop.is_Symbol:
                yield node + ', ' + str(sop)
            else:
                if sop.is_Function:
                    time.sleep(10)
                    yield node + ', ' + '|'.join(str(sop).replace(' ', '').replace('~', '!').split('|'))
                else:
                    (yield node + ', 1') if sop.canonical else (yield node + ', 0')

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

        if 'space' not in dir(self):
            combinations = [bin(i).split('b')[1] if len(bin(i).split('b')[1]) == len(self.graph.index) else
                            ''.join(str_gen(len(self.graph.index) - len(bin(i).split('b')[1]))) + bin(i).split('b')[1]
                            for i in range(0, 2 ** len(self.graph.index))]
            self.space = [{self.graph.index[i]: int(c[i]) for i in range(0, len(c))} for c in combinations]

        return self.space

    def get_support(self, node=None, hex_flag=False):
        """
        DESCRIPTION:
        A method to obtain the support of a given node. If not, it returns the support of all the nodes.
        :param node: [string] name of the node whose support must be obtained.
        :param hex_flag: [boolean] flag to determine whether the support is to be given in hexadecimal or not.
        :return: [dict] all the terms of support selected by node.
        """
        # Auxiliary functions
        def converter(number):
            if hex_flag:
                return (lambda x: hex(int(x, 2)))(number)
            else:
                return (lambda x: x)(number)

        if node is not None:
            return {node: list(map(converter, self.maps.loc[node, :][lambda x: x].index.values.tolist()))}
        else:
            return {node: list(map(converter, self.maps.loc[node, :][lambda x: x].index.values.tolist()))
                    for node in self.maps.index}

    def set_maps(self):
        """
        DESCRIPTION:
        A method to calculate the response of the function in the whole space of the boolean expressions of the network.
        :return: [pandas DataFrame] the maps for all the nodes.
        """
        # Auxiliary functions
        def get_symbol(source, letter):
            """
            DESCRIPTION:
            It is here where the role is introduced.
            """
            symbol = '!'
            if self.roles_set is not None:
                canalizing_value = list(filter(lambda x: x[0] == source and x[1] == letter, self.roles_set))[0][2]
                symbol = '!' if canalizing_value == 1 else ''
            return symbol

        # Build the expression
        network = []
        for i in range(0, len(self.network)):
            # Set the structure of an INPUT node
            if self.network[i] == 'INPUT':
                network.append(f'!(!{self.graph.index[i]})')
                continue
            # Iterate through the layers of a node
            load = ''
            steps = range(len(self.network[i]) - 1, -1, -1)
            m_step = min(steps)
            first_letter = list(self.network[i][m_step])[0]
            for j in steps:
                layer = self.network[i][j]
                # Calculate the whole load
                if len(layer) == 1:
                    level = get_symbol(source=self.graph.index[i], letter=layer[0]) + layer[0]
                else:
                    level = get_symbol(source=self.graph.index[i], letter=layer[0]) + layer[0]
                    for letter in layer[1::]:
                        level = get_symbol(source=self.graph.index[i], letter=letter) + letter + '&' + level
                # Add inner load
                level = level + '&' + load if load != '' else level
                # Add the header
                if j != m_step or j == m_step and first_letter in self.graph.loc[self.graph.index[i], 'activators']:
                    load = '!(' + level + ')'
                else:
                    load = level
            network.append(load)
        # Build the map
        values = self.eval_expression(self.get_space(), network)
        values = [values[int(i*len(values)/len(self.graph.index)):int((i+1)*len(values)/len(self.graph.index))] for i in range(0, len(self.graph.index))]
        # Change and send the map
        kmaps = pd.DataFrame(data=values, index=self.graph.index, columns=list(map(
            lambda x: bin(x)[2::].zfill(len(self.graph.index)), range(0, len(self.get_space())))))
        self.maps = kmaps
        return kmaps


def simplify(symbols, minterms):
    """
    DESCRIPTION:
    A function that, given a set of minterm in binary form, simplifies the whole set into the minimal sum of products.
    :param symbols: [list] variables (strings) representing the nodes.
    :param minterms: [list] the whole group (ints) of minterms to be simplified. The numbers represent the binary order.
    :return: [Symbol] products of the sum of products in boolnet format.
    """
    variables = ' '.join(symbols)
    return SOPform(variables=sympy.symbols(variables), minterms=minterms, dontcares=[])






