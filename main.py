
"""
DESCRIPTION: 
In this script we test all the functions developed. 
"""

import numpy as np 
import pandas as pd
from utils.ncbf import ncbfCalc 

def main():
    # Nodes: X, A, B, C.
    # Build the input object.
    index = np.array(['X', 'A', 'B', 'C'])
    columns = np.array(['activators', 'inhibitors'])
    data = np.array([[['A', 'B', 'C', 'D'], ['E', 'F', 'G']], [['E', 'X', 'G'], ['A', 'B', 'C']], [['A', 'X', 'G'], ['F', 'B', 'C', 'E']], [['G', 'D', 'X'], ['A', 'F', 'B']]])
    data = pd.DataFrame(data=data, index=index, columns=columns)

    # We get another DataFrame with all the possible ncbf
    paths = ncbfCalc(data=data)

main()