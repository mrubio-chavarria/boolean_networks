
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
    index = ['I', 'S', 'T', 'Z', 'D']
    columns = ['activators', 'inhibitors']
    content = [[[''], ['']], [['I'], ['S', 'T']], [[''], ['S', 'Z']], [['S', 'Z'], ['D']], [['D'], ['S', 'Z']]]
    initial_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Generate the graph object
    attractors = ['10010', '11101', '01101']  # Introduce each one in alphabetical order
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
                                    ['D', 'Z', 0, 0],
                                    ['D', 'D', 0, 0]]],
                     'structures': [['INPUT', ['I', 'S', 'T'], ['SZ'], ['D', 'SZ'], ['S', 'ZD']]]
                     }
    imposed_roles_sets = filter_kernel['roles_sets']
    graph = Graph(initial_data=initial_data, attractors=attractors, filter_kernel=filter_kernel,
                  imposed_roles_sets=imposed_roles_sets)

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

