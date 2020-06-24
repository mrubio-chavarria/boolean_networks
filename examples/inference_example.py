#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
INDICATIONS:

    In this script it is reproduced and automatized the inference of the epithelial-to-mesenchymal transition (EMT). The
biological information has  been taken from:

    - Joo, J.I., Zhou, J.X., Huang, S. et al. Determining Relative Dynamic Stability of Cell States Using Boolean
      Network Model. Sci Rep 8, 12077 (2018). https://doi.org/10.1038/s41598-018-30544-0

    - Modeling EMT decision circuit Mingyang Lu, Mohit Kumar Jolly, Herbert Levine, José N. Onuchic, Eshel Ben-Jacob
      Proceedings of the National Academy of Sciences Nov 2013, 110 (45) 18144-18149; DOI: 10.1073/pnas.1318192110

"""

from graphs.graph import Graph
import pandas as pd


def main():
    # ------------------------------------------------------------------------------------------------------------------
    # EMT
    # ------------------------------------------------------------------------------------------------------------------
    # Introduce the data
    index = ['A', 'B', 'C', 'D', 'E']
    columns = ['activators', 'inhibitors']
    content = [[[''], ['']], [['A'], ['B', 'C']], [[''], ['B', 'D']], [['B', 'D'], ['E']], [[''], ['B', 'D']]]
    initial_data = pd.DataFrame(data=content, index=index, columns=columns)
    attractors = ['00101', '11011', '11010']  # Introduce each one in alphabetical order

    # Inference parameters
    simulations = 20
    variants_limit = None
    max_local_iterations = 50
    max_global_iterations = 50
    filter_kernel = {'roles_sets': [[['A', 'A', 1, 1],
                                     ['B', 'A', 0, 0],
                                     ['B', 'B', 0, 1],
                                     ['B', 'C', 0, 1],
                                     ['C', 'B', 1, 0],
                                     ['C', 'D', 1, 0],
                                     ['D', 'E', 0, 1],
                                     ['D', 'B', 0, 0],
                                     ['D', 'D', 0, 0],
                                     ['E', 'B', 0, 1],
                                     ['E', 'D', 0, 1]]],
                     'structures': [[['BD'], 'INPUT', ['A', 'BC'], ['BD'], ['E', 'BD']]]}
    imposed_roles_sets = None

    # Model inference
    graph = Graph(initial_data=initial_data, attractors=attractors, filter_kernel=None,
                  imposed_roles_sets=imposed_roles_sets, simulations=simulations, variants_limit=variants_limit,
                  max_global_iterations=max_global_iterations, max_local_iterations=max_local_iterations)

    # Print results in file
    file = open('results.txt', 'r+')
    file.truncate(0)
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nTotal results\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['total_results']]
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nNon-cyclic results\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['non_cyclic_results']]
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nSame number of attractors\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['same_number']]
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nOne attractor of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['one_of_all']]
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nTwo attractors of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['two_of_all']]
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nThree attractors of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['three_of_all']]
    file.close()


if __name__ == '__main__':
    main()



