
"""
DESCRIPTION: 
In this script we test all the functions developed. 
"""

import numpy as np 
import pandas as pd
from utils.ncbf import ncbfCalc, networksCalc, netValidator

def main():

    # Build the input object
    index = ['A', 'B', 'C', 'D']
    columns = ['activators', 'inhibitors']
    data = [[['D'], ['']], [['A'], ['B', 'C']], [['B', 'C'], ['A']], [['B'], ['C']]]
    data = pd.DataFrame(data=data, index=index, columns=columns)

    # Get another DataFrame with all the possible ncbf. TO BE IMPROVED.
    paths = ncbfCalc(data=data)

    # Get all possible networks given the set of paths
    networks = list(networksCalc(paths))

    # Return the set of networks which meet the condition
    attractors = [[0, 0, 1, 0], [1, 1, 0, 1]]
    networks = netValidator(networks, data, attractors)

if __name__ == '__main__':
    main()