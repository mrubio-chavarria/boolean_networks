
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
    # Example of the article of Huan 2018
    # ------------------------------------------------------------------------------------------------------------------
    # Introduce the data
    index = ['I', 'S', 'T', 'Z', 'D']
    columns = ['activators', 'inhibitors']
    content = [[[''], ['']], [['I'], ['S', 'T']], [[''], ['S', 'Z']], [['S', 'Z'], ['D']], [['D'], ['S', 'Z']]]
    initial_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Generate the graph object
    graph = Graph(initial_data=initial_data)

    # Get another DataFrame with all the possible ncbf.
    tags = ['activators', 'inhibitors']
    initial_paths = ncbfCalc(data=initial_data, tags=tags)

    # Get all possible networks given the set of paths
    initial_networks = list(networksCalc(initial_paths))

    # Return the set of networks which meet the condition
    attractors = ['00101', '11011', '11010']
    final_networks = netValidator(initial_networks=initial_networks, initial_graph=initial_data, tags=tags,
                                  attractors=attractors)

    # Filter by networks with coherent projection
    networks = netFilter(final_networks)


if __name__ == '__main__':
    main()

