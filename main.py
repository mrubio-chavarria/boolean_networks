
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
    attractors = ['0010', '1101']  # Introduce each one in alphabetical order
    # Generate the graph object
    simulations = 2
    variants_limit = 3
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
    # imposed_roles_sets = filter_kernel['roles_sets']
    graph = Graph(initial_data=initial_data, attractors=attractors, filter_kernel=None,
                  imposed_roles_sets=None, simulations=simulations, variants_limit=variants_limit,
                  max_global_iterations=max_global_iterations, max_local_iterations=max_local_iterations)
    """
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
    """

if __name__ == '__main__':
    main()

