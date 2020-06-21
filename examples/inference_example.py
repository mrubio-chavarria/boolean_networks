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
    index = ['A', 'B', 'C', 'D']
    columns = ['activators', 'inhibitors']
    content = [[['D'], ['']], [['A'], ['']], [[''], ['A']], [['B'], ['C']]]
    initial_data = pd.DataFrame(data=content, index=index, columns=columns)
    attractors = ['0010', '1101']  # Introduce each one in alphabetical order

    # Inference parameters
    simulations = 1
    variants_limit = 1
    max_local_iterations = 50
    max_global_iterations = 100
    filter_kernel = {'roles_sets': [[['I', 'I', 1, 1],
                                    ['S', 'I', 0, 0],
                                    ['S', 'S', 0, 1],
                                    ['S', 'T', 0, 1],
                                    ['T', 'S', 1, 0],
                                    ['T', 'Z', 1, 0],
                                    ['Z', 'D', 0, 1],
                                    ['Z', 'S', 0, 0],
                                    ['Z', 'Z', 0, 0],
                                    ['D', 'S', 0, 1],
                                    ['D', 'Z', 0, 1]]],
                     'structures': [[['SZ'], 'INPUT', ['I', 'ST'], ['SZ'], ['D', 'SZ']]]
                     }
    imposed_roles_sets = None

    # Model inference
    graph = Graph(initial_data=initial_data, attractors=attractors, filter_kernel=None,
                  imposed_roles_sets=imposed_roles_sets, simulations=simulations, variants_limit=variants_limit,
                  max_global_iterations=max_global_iterations, max_local_iterations=max_local_iterations)


if __name__ == '__main__':
    main()



