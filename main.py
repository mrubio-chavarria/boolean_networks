
"""
DESCRIPTION: 
In this script we test all the functions developed. 
"""

from utils.graph import Graph
import pandas as pd
from utils.ncbf import ncbfCalc, conflict_ncbfCalc, networksCalc, netValidator, netFilter
import ast


def main():
    # ------------------------------------------------------------------------------------------------------------------
    # EMT
    # ------------------------------------------------------------------------------------------------------------------
    # Example of the article of Huang2018
    # ------------------------------------------------------------------------------------------------------------------
    # Introduce the data
    index = ['A', 'B', 'C', 'D']
    columns = ['activators', 'inhibitors']
    content = [[['D'], ['']], [['A'], ['']], [[''], ['A']], [['B'], ['C']]]
    initial_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Generate the graph object
    simulations = 20
    variants_limit = None
    max_local_iterations = 100
    max_global_iterations = 200
    attractors = ['1101', '0010']  # Introduce each one in alphabetical order
    filter_kernel = {'roles_sets': [[['I', 'I', 1, 1],
                                    ['S', 'I', 1, 1],
                                    ['S', 'S', 0, 0],
                                    ['S', 'T', 0, 1],
                                    ['T', 'S', 1, 0],
                                    ['T', 'Z', 1, 0],
                                    ['Z', 'D', 0, 1],
                                    ['Z', 'S', 0, 0],
                                    ['Z', 'Z', 0, 0],
                                    ['D', 'S', 0, 1],
                                    ['D', 'Z', 0, 1]]],
                     'structures': [['INPUT', ['I', 'S', 'T'], ['SZ'], ['D', 'SZ'], ['SZ']]]
                     }
    filter_kernel = {'roles_sets': [[['A', 'D', 1, 1],
                                     ['B', 'A', 1, 1],
                                     ['C', 'A', 1, 0],
                                     ['D', 'B', 1, 1],
                                     ['D', 'C', 1, 0]]],
                     }
    # imposed_roles_sets = filter_kernel['roles_sets']
    graph = Graph(initial_data=initial_data, attractors=attractors, filter_kernel=filter_kernel,
                  imposed_roles_sets=None, simulations=simulations, variants_limit=None,
                  max_global_iterations=max_global_iterations, max_local_iterations=max_local_iterations)

    # Get another DataFrame with all the possible ncbf.
    tags = ['activators', 'inhibitors']
    initial_paths = ncbfCalc(data=initial_data, tags=tags)

    # Get all possible networks given the set of paths
    initial_networks = list(networksCalc(initial_paths))

    # Return the set of networks which meet the condition
    final_networks = netValidator(initial_networks=initial_networks, initial_graph=initial_data, tags=tags,
                                  attractors=attractors)

    # Filter by networks with coherent projection
    networks = netFilter(final_networks)


if __name__ == '__main__':
    main()

