#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import sys
from uuid import uuid1


class Result:
    """
    DESCRIPTION:
    A class built to store the results of the simulation.
    """

    def __init__(self, network, pathways, maps_set, conflicts, simulation, iterations, expressions, attractors, variant,
                 priority_matrix, roles_set, accepted=False):
        """
        DESCRIPTION:
        Constructor to the Result class.
        :param network: [list] the structured network employed.
        :param pathways: [list] objects pathway, the result of the simulation.
        :param conflicts: [list] objects conflict, the record of the simulation.
        :param maps_set: [KMap] the set of maps employed during the simulation.
        :param simulation: [int] identifier of the simulation.
        :param iterations: [int] number of the iterations performed.
        :param expressions: [list] final expressions of the net in boolnet format.
        :param attractors: [dict] all the attractors of the network, steady and cyclic.
        :param variant: [Variant] the modification of the graph with the roles that produced this result.
        :param priority_matrix: [pandas DataFrame] a matrix to represent to priorities taken.
        :param roles_set: [list] group of pairs canalized-canalizing values for expressions.
        :param accepted: [boolean] sign to indicate whether the pathways has been included or not.
        """
        self.network = network
        self.pathways = None
        self.set_pathways(pathways)
        self.maps_set = maps_set
        self.conflicts = conflicts
        self.code = None
        self.code = self.get_code()
        self.iterations = iterations
        self.simulation = simulation
        self.expressions = expressions
        self.attractors = attractors
        self.accepted = accepted
        self.variant = variant
        self.priority_matrix = priority_matrix
        self.roles_set = roles_set

    def get_code(self):
        """
        DESCRIPTION:
        A method to built a unique identifier of the type to which the result belongs.
        :return: [string] code of the kind of pathway
        """
        if self.code is None:
            codes = ['$']
            [codes.append(f'{str(hex(int("".join(pathway.region_of_interest), 2)))}:{pathway.consequent}|')
             for pathway in self.pathways]
            code = ''.join(codes) + '$'
        else:
            code = self.code
        return code

    def get_serialized_data(self):
        """
        DESCRIPTION:
        A method to return the serialized value of the result. A readable representation of the result to be introduced
        in a file.
        :return: [dictionary] representation of the result in JSON format.
        """
        data = {'network': str(self.network),
                'pathways': '\n                  '.join([str(pathway) for pathway in self.pathways]),
                'conflicts': '\n                   '.join([str(conflict) for conflict in self.conflicts]),
                'attractors': '\n                    ' +
                              'Steady: ' + ' '.join(self.attractors['steady']) +
                              '\n                    ' +
                              'Cyclic: ' + ' '.join([str(att) for att in self.attractors['cyclic']]),
                'accepted': str(self.accepted),
                'matrix': '\n        ' + '\n        '.join(str(self.priority_matrix).split('\n')),
                'roles': '\n        ' + str(self.variant.roles),
                'expressions': '\n        ' + '\n        '.join(self.expressions),
                'final_map': '\n        ' + '\n        '.join(str(self.maps_set.maps).split('\n'))}
        text = f"""
        * Network:
        ID: {uuid1()}
        Structure: {data['network']}
        Pathways: {data['pathways']}
        Conflicts: {data['conflicts']}
        Attractors: {data['attractors']}
        Accepted: {data['accepted']}
        Priority matrix: {data['matrix']}
        Roles set: {data['roles']}
        Expressions: {data['expressions']}
        Final map: {data['final_map']}
        """
        return text

    def set_pathways(self, pathways):
        """
        DESCRIPTION:
        A method to get the pathways always ordered in the same manner.
        :param pathways: [list] objects pathway, the result of the simulation.
        """
        pathways.sort(key=lambda x: x.consequent)
        self.pathways = pathways



