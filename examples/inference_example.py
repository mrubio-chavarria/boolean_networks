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
    simulations = 1
    variants_limit = None
    max_local_iterations = 50
    max_global_iterations = 50
    roles_sets = [[['A', 'A', 1, 1],
                   ['B', 'A', 0, 0],
                   ['B', 'B', 0, 1],
                   ['B', 'C', 0, 1],
                   ['C', 'B', 1, 0],
                   ['C', 'D', 1, 0],
                   ['D', 'E', 0, 1],
                   ['D', 'B', 0, 0],
                   ['D', 'D', 0, 0],
                   ['E', 'B', 0, 1],
                   ['E', 'D', 0, 1]]]
    structures = [['INPUT', ['ABC'], ['BD'], ['D', 'B', 'E'], ['B', 'D']],
                  ['INPUT', ['AB', 'C'], ['BD'], ['B', 'DE'], ['BD']],
                  ['INPUT', ['ABC'], ['BD'], ['D', 'B', 'E'], ['D', 'B']],
                  ['INPUT', ['A', 'C', 'B'], ['B', 'D'], ['B', 'D', 'E'], ['BD']]]
    filter_kernel = {'roles_sets': None,
                     'structures': structures}
    imposed_roles_sets = None

    # Model inference
    graph = Graph(initial_data=initial_data, attractors=attractors, filter_kernel=filter_kernel,
                  imposed_roles_sets=imposed_roles_sets, simulations=simulations, variants_limit=variants_limit,
                  max_global_iterations=max_global_iterations, max_local_iterations=max_local_iterations)

    # Print results in files
    file = open('results_non_cyclic.txt', 'w')
    file.truncate(0)
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nNon-cyclic results\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['non_cyclic_results']]
    file.close()
    file = open('results_same_number.txt', 'w')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nSame number of attractors\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['same_number']]
    file.close()
    file = open('results_one.txt', 'w')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nOne attractor of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['one_of_all']]
    file.close()
    file = open('results_same_one.txt', 'w')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nSame number of attractors and one attractor of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['same_number_with_at_least_one_of_all']]
    file.close()
    file = open('results_two.txt', 'w')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nTwo attractors of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['two_of_all']]
    file.close()
    file = open('results_same_two.txt', 'w')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nSame number of attractors and two attractors of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['same_number_with_at_least_two_of_all']]
    file.close()
    file = open('results_three.txt', 'w')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nThree attractors of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['three_of_all']]
    file.close()
    file = open('results_same_three.txt', 'w')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    file.write('\nSame number and three attractors of all\n')
    file.write('\n--------------------------------------------------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in graph.results['same_number_and_three_of_all']]
    file.close()


if __name__ == '__main__':
    main()



