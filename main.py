
"""
DESCRIPTION: 
In this script we test all the functions developed. 
"""

import numpy as np 
import pandas as pd
from utils.ncbf import ncbfCalc, networksCalc

def main():

    # Build the input object.
    index = np.array(['I', 'A', 'B', 'C', 'D'])
    columns = np.array(['activators', 'inhibitors'])
    data = np.array([[[''], ['']], [['I'], ['C']], [['A'], ['D']], [['B'], ['']], [['B', 'C'], ['A']]])
    data = pd.DataFrame(data=data, index=index, columns=columns)

    # We get another DataFrame with all the possible ncbf.
    paths = ncbfCalc(data=data)

    # We get all possible networks given the set of paths.
    networks = list(networksCalc(paths))
    print(networks)

main()