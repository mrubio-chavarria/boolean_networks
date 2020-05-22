
"""
DESCRIPTION: 
In this script we test all the functions developed. 
"""

import numpy as np 
import pandas as pd
from utils.ncbf import ncbfCalc, conflict_ncbfCalc, networksCalc, netValidator

def main():

    # Original graph
    index = ['A', 'B', 'C', 'D']
    columns = ['activators', 'inhibitors']
    content = [[['D', 'F'], ['']], [['A'], ['E']], [['E'], ['A']], [['B'], ['C']]]
    original_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Conflicts graph
    index = ['E', 'F']
    tag = 'variables'
    columns = [tag]
    content = [[['A', 'E']], [['A', 'B', 'C']]]
    conflicts_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Get another DataFrame with all the possible ncbf.
    tags = ['activators', 'inhibitors']
    original_paths = ncbfCalc(data=original_data, tags=tags)
    conflicts_paths, conflicts_graphs = conflict_ncbfCalc(variables=conflicts_data, tag=tag)

    # Get all possible networks given the set of paths
    original_networks = list(networksCalc(original_paths))
    sets_conflicts_networks = []
    for set_paths in conflicts_paths:
        sets_conflicts_networks.append(list(networksCalc(set_paths)))

    # Return the set of networks which meet the condition
    attractors = ['1101', '0010']
    final_networks = []
    num_sets_conflicts_networks = len(sets_conflicts_networks)
    for i in range(0, num_sets_conflicts_networks):
        final_networks.extend(
            netValidator(original_networks, original_data, attractors, sets_conflicts_networks[i], conflicts_graphs[i], tags)
        )

    # Save the obtained set of networks
    file = open('networks.txt', 'w')
    file.write(str(final_networks))
    file.close()
    print()

if __name__ == '__main__':
    main()

