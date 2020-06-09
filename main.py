
"""
DESCRIPTION: 
In this script we test all the functions developed. 
"""

import numpy as np
import pandas as pd
from utils.ncbf import ncbfCalc, conflict_ncbfCalc, networksCalc, netValidator, netFilter
import ast


def main():

    # Original graph
    index = ['A', 'B', 'C', 'D']
    columns = ['activators', 'inhibitors']
    content = [[['D', 'F'], ['']], [['A'], ['E']], [['E'], ['A']], [['B'], ['C']]]
    original_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Fixed conflicts graphs.
    # Be noticed that here they are the graph with the variables and the one without any conflict introduced
    index = ['E', 'F']
    columns = ['expressions', 'activated', 'inhibited']
    content = [[['B', 'C'], ['C'], ['B']], [['A', 'B', 'C'], ['A'], ['']]]
    fixed_conflicts_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Introduce the example of Huang
    index = ['I', 'S', 'T', 'Z', 'D']
    columns = ['activators', 'inhibitors']
    content = [[[''], ['']], [['I'], ['S', 'T']], [[''], ['S', 'Z']], [['S', 'Z'], ['D']], [['D'], ['S', 'Z']]]
    initial_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Unfixed conflicts graphs
    index = ['E', 'F']
    tag = 'variables'
    columns = [tag]
    content = [[['A', 'E']], [['A', 'D', 'E', 'F']]]
    unfixed_conflicts_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Get another DataFrame with all the possible ncbf.
    tags = ['activators', 'inhibitors']
    initial_paths = ncbfCalc(data=initial_data, tags=tags)
    original_paths = ncbfCalc(data=original_data, tags=tags)
    unfixed_conflicts_paths, unfixed_conflicts_graphs = conflict_ncbfCalc(variables=unfixed_conflicts_data, tag=tag)

    # Get all possible networks given the set of paths
    original_networks = list(networksCalc(original_paths))
    initial_networks = list(networksCalc(initial_paths))
    unfixed_sets_conflicts_networks = []
    for set_paths in unfixed_conflicts_paths:
        unfixed_sets_conflicts_networks.append(list(networksCalc(set_paths)))

    # Return the set of networks which meet the condition
    attractors = ['00101', '11011', '11010']
    final_networks = netValidator(initial_networks, initial_data, original_networks, original_data, attractors,
                                  unfixed_sets_conflicts_networks, unfixed_conflicts_graphs, fixed_conflicts_data, tags)

    # Filter by networks with coherent projection
    networks = netFilter(final_networks)


if __name__ == '__main__':
    main()

