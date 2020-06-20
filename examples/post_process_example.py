#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DESCRIPTION:
In this script it is automatized the inference of the example network described in the dissertation.
"""

from base.ncbf import ncbfCalc, conflict_ncbfCalc, networksCalc, post_process
import pandas as pd


def main():
    # ------------------------------------------------------------------------------------------------------------------
    # Extension of the example network
    # ------------------------------------------------------------------------------------------------------------------
    # Functional extension with the conflicts of the network
    # ------------------------------------------------------------------------------------------------------------------
    # Original graph
    index = ['A', 'B', 'C', 'D']
    columns = ['activators', 'inhibitors']
    content = [[['D', 'F'], ['']], [['A'], ['E']], [['E'], ['A']], [['B'], ['C']]]
    original_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Unfixed conflicts graphs
    index = ['E', 'F']
    tag = 'variables'
    columns = [tag]
    content = [[['A', 'E']], [['A', 'D', 'E', 'F']]]
    unfixed_conflicts_data = pd.DataFrame(data=content, index=index, columns=columns)

    # Get another DataFrame with all the possible ncbf.
    tags = ['activators', 'inhibitors']
    original_paths = ncbfCalc(data=original_data, tags=tags)
    unfixed_conflicts_paths, unfixed_conflicts_graphs = conflict_ncbfCalc(variables=unfixed_conflicts_data, tag=tag)

    # Get all possible networks given the set of paths
    original_networks = list(networksCalc(original_paths))
    unfixed_sets_conflicts_networks = []
    for set_paths in unfixed_conflicts_paths:
        unfixed_sets_conflicts_networks.append(list(networksCalc(set_paths)))

    # Return the set of networks which meet the condition
    attractors = ['00101', '11011', '11010']
    results = post_process(unfixed_sets_conflicts_networks=unfixed_sets_conflicts_networks,
                           unfixed_conflicts_graphs=unfixed_conflicts_graphs,
                           original_networks=original_networks,
                           original_graph=original_data,
                           attractors=attractors)


if __name__ == '__main__':
    main()

