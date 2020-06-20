
"""
DESCRIPTION: 
In this script we test all the functions developed. 
"""

from graphs.graph import Graph
import pandas as pd


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
    # ------------------------------------------------------------------------------------------------------------------
    # Generate the graph object
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
    # imposed_roles_sets = filter_kernel['roles_sets']
    graph = Graph(initial_data=initial_data, attractors=attractors, filter_kernel=None,
                  imposed_roles_sets=None, simulations=simulations, variants_limit=variants_limit,
                  max_global_iterations=max_global_iterations, max_local_iterations=max_local_iterations)
    """
    # ------------------------------------------------------------------------------------------------------------------
    # Functional extension with the conflicts of the network
    # ------------------------------------------------------------------------------------------------------------------
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

    # Unfixed conflicts graphs
    index = ['E', 'F']
    tag = 'variables'
    columns = [tag]
    content = [[['A', 'E']], [['A', 'D', 'E', 'F']]]
    unfixed_conflicts_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Introduce the example of Huang
    index = ['I', 'S', 'T', 'Z', 'D']
    columns = ['activators', 'inhibitors']
    content = [[[''], ['']], [['I'], ['S', 'T']], [[''], ['S', 'Z']], [['S', 'Z'], ['D']], [['D'], ['S', 'Z']]]
    initial_data = pd.DataFrame(data=content, index=index, columns=columns)

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
    results = third_validation(graph=initial_data, initial_networks=initial_networks, attractors=attractors)
    final_networks = netValidator(initial_networks=initial_networks, initial_graph=initial_data,
                                  original_networks=original_networks, original_graph=original_data,
                                  attractors=attractors, unfixed_sets_conflicts_networks=unfixed_sets_conflicts_networks,
                                  unfixed_conflicts_graphs=unfixed_conflicts_graphs, tags=tags)

    # Filter by networks with coherent projection
    networks = netFilter(final_networks)
    """

if __name__ == '__main__':
    main()

