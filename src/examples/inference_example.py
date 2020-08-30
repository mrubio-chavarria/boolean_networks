#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from graphs.graph import Graph, Manager
import pandas as pd
import os
import sys


def main():

    # Introduce the data
    index = ['A', 'B', 'L', 'M']
    columns = ['activators', 'inhibitors']
    content = [[['L', 'B'], ['A', 'B']], [['M'], ['B']], [[''], ['']], [['A'], ['M']]]
    initial_data = pd.DataFrame(data=content, index=index, columns=columns)
    attractors = ['0000', '1111']  # Introduce each one in alphabetical order

    # Inference parameters
    simulations = 200
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
    structures = [[['BL', 'A'], ['M', 'B'], ['M', 'A'], 'INPUT']]
    filter_kernel = {'roles_sets': None,
                     'structures': None}
    imposed_roles_sets = [[['A', 'A', 0, 1, '01'], ['A', 'B', 1, 0, '10'], ['A', 'L', 0, 0, '00'], ['B', 'B', 0, 1, '01'], ['B', 'M', 0, 0, '00'], ['L', 'L', 1, 1, '11'], ['M', 'A', 1, 1, '11'], ['M', 'M', 1, 0, '10']]]

    # Model inference
    manager = Manager(
        initial_data=initial_data, attractors=attractors,
        filter_kernel=filter_kernel, imposed_roles_sets=None,
        simulations=simulations, variants_limit=variants_limit,
        max_global_iterations=max_global_iterations,
        max_local_iterations=max_local_iterations
    )

    # Print results in files
    file = open('results_cyclic.txt', 'w')
    file.truncate(0)
    file.write('\n------------------ --------------------------------------\n')
    file.write('\nCyclic results\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['cyclic_results']]
    file.close()
    file = open('results_non_cyclic.txt', 'w')
    file.truncate(0)
    file.write('\n---------------------------------------------------------\n')
    file.write('\nNon-cyclic results\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['non_cyclic_results']]
    file.close()
    file = open('results_same_number.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nSame number of attractors\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['same_number']]
    file.close()
    file = open('results_one.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nOne attractor of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['one_of_all']]
    file.close()
    file = open('results_same_one.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nSame number of attractors and one attractor of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['same_number_with_at_least_one_of_all']]
    file.close()
    file = open('results_two.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nTwo attractors of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['two_of_all']]
    file.close()
    file = open('results_same_two.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nSame number of attractors and two attractors of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['same_number_with_at_least_two_of_all']]
    file.close()
    file = open('results_three.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nThree attractors of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['three_of_all']]
    file.close()
    file = open('results_same_three.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nSame number and three attractors of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['same_number_and_three_of_all']]
    file = open('results_four.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nFour attractors of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['four_of_all']]
    file.close()
    file = open('results_same_four.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nSame number and four attractors of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['same_number_and_four_of_all']]
    file.close()
    file = open('results_same_four.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nSame number and four attractors of all\n')
    file.write('\n---------------------------------------------------------\n')
    [file.write('\n' + result.get_serialized_data()) for result in manager.results['same_number_and_four_of_all']]
    file.close()
    # Print the analysis file
    file = open('analysis.txt', 'w')
    file.write('\n---------------------------------------------------------\n')
    file.write('\nAnalysis performed\n')
    file.write('\n---------------------------------------------------------\n')
    file.write(manager.get_serialized_analysis())
    file.close()

if __name__ == '__main__':
    main()
