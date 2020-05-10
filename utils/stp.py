

"""
INDICATIONS:
In this file, they are given several utilities to construct boolean networks based on the Semi-Tensor product. It is all
based in the structure stated in Murrugarra 2013.
"""

import numpy as np
from numpy import matlib as mb
from random import choice

def lm(columns, n):
    """
    DESCRIPTION:
    A function devised to build logic matrices based on delta notation.
    :param columns: [numpy array] array with the indexes of identity columns.
    :param n: [int] dimension of the nxn identity matrix.
    """
    columns = columns - np.ones(columns.shape)
    i = np.eye(n)
    e = i[:, columns.astype(int).tolist()]
    return e




def dumM(k=2):
    """
    DESCRIPTION:
    A function to generate the dummy operand with k dimensions. The standard dummy operand is for k = 2. Adapted from
    the one found in Analysis and Control of Boolean Networks.
    :param k: [int] dimension to create the dummy operand.
    :return: [numpy array] dummy operand.
    """

    a = np.linspace(1, k, k)
    return lm(mb.repmat(a, 1, k)[0, :], k)


def stp(a, b):
    """
    DESCRIPTION:
    A function to compute th left semi-tensor product of two matrices.
    :param a: [numpy array] first matrix to multiply.
    :param b: [numpy array] second matrix to multiply.
    :return: [numpy array] the semi-tensor product of a and b.
    """

    # Obtain the dimensions of boths matrices
    n = a.shape[1]
    p = b.shape[0]
    c = np.lcm(n, p)

    # Obtain the identity matrices
    ia = np.eye(int(c/n))
    ib = np.eye(int(c/p))

    # Compute the semi-tensor product of both matrices
    sp = np.dot(np.kron(a, ia), np.kron(b, ib))

    return sp

def swapM(m, n):
    """
    DESCRIPTION:
    A function that creates a swap matrix provided the dimensions for the indexes. Adapted from the one found in
    Analysis and Control of Boolean Networks.
    :param m: [int] first dimension.
    :param n: [int] second dimension.
    :return: [numpy array] the swap matrix.
    """

    # Parameters
    d = m*n
    w = np.zeros((d, d), dtype=int)

    # Compute the swap matrix
    for k in range(0, d):
        j = (k+1) % n
        if j == 0:
            j = n
        i = (k-j+1)/n + 1
        e = int((j-1)*m+i)
        w[e-1, k] = 1

    return w


def phiGen(s):
    """
    DESCRIPTION:
    A function to generate the complementary factor in the formula to generate the matrix of a boolean network according
    with the algebraic form.
    :param s: [int] the subindex of phi function.
    :return: [numpy array] the result of the phi function.
    """

    # Parameters
    mr = np.array([[1, 0], [0, 0], [0, 0], [0, 1]], dtype=int)  # Power reducing matrix
    base = lambda n, j: np.kron(np.eye(2 ** (j-1)), stp(np.kron(np.eye(2), swapM(2, 2 ** (n-j))), mr))
    load = base(s, s)

    # Compute phi
    for i in range(s-1, 0, -1):
        load = stp(base(s, i), load)
    phi = load

    return phi


def convertSP(vars):
    if len(vars) == 1:
        return vars[0]
    else:
        chunk = '*'.join(vars)
        return chunk


def convertKP(vars, logic):
    return f'(leye({logic})+{vars})'


def swapvm(vars, logic):
    value1 = convertSP(vars[1::])
    value2 = convertKP(value1, logic)
    return [value2, vars[0]]


def PR_Swap(vars, f, logic):
    """
    DESCRIPTION:
    MR: power reduction matrix.
    """
    for i in range(1, max(f)+1):
        v = [idx[0] + 1 for idx in enumerate(f) if idx[1] == i]
        j = 2
        while j <= len(v):
            chunk = f[v[j - 2]: v[j - 1] - 1]
            if v[j-1]-v[j-2] == 1:
                vars[v[j-2]-1] = 'MR'
                f[v[j-2]-1] = 0
                v = [idx[0] + 1 for idx in enumerate(f) if idx[1] == i]
                continue
            elif not any(chunk):
                vars[v[j-2]] = convertKP(convertSP(vars[v[j-2]:v[j-1]-1]), logic)
                f[v[j-2]-1: v[j-2] + 1] = [0 for item in f[v[j-2]-1: v[j-2] + 1]]
                vars[v[j - 2]] = 'MR'
                k = range(v[j - 2] + 1, v[j-1] - 1, 1)
                positions = [u for u in range(0, len(f)) if u not in [c - 1 for c in k]]
                f = [f[pos] for pos in positions]
                vars = [vars[pos] for pos in positions]
                v = [idx[0] + 1 for idx in enumerate(f) if idx[1] == i]
                continue
            j = j + 1
    return vars, f


def moveVarBack(vars, f, logic):
    ln = len(f)
    v = [idx[0] + 1 for idx in enumerate(f) if idx[1] > 0]
    l = len(v)
    for i in range(l, 0, -1):
        if v[i-1] == ln:
            continue
        chunk = [idx[0] + 1 for idx in enumerate(f[v[i-1]:ln-(l-i)]) if idx[1] == 0]
        if any(chunk):
            vars[v[i-1]-1: v[i-1] + 1] = swapvm(vars[v[i-1]-1: ln - (l - i)], logic)
            f[v[i-1]] = f[v[i-1]-1]
            f[v[i-1]-1] = 0
            if v[i-1] + 2 - (ln - (l - i)) == 0:
                k = [v[i-1] + 2]
            else:
                k = range(v[i-1] + 2, ln - (l - i) + 1, 1)
            if len(k) != 0:
                positions = [u for u in range(0, len(f)) if u not in [c-1 for c in k]]
                f = [f[pos] for pos in positions]
                vars = [vars[pos] for pos in positions]
            ln = len(f)
    return vars, f

def swapVarsByOrder(vars, f, logic):

    k = logic
    v = [idx[0] + 1 for idx in enumerate(f) if idx[1] == 0]
    vars1 = [vars[i-1] for i in v]
    f1 = [f[i-1] for i in v]
    positions = [u for u in range(0, len(f)) if u not in [c - 1 for c in v]]
    vars = [vars[pos] for pos in positions]
    f = [f[pos] for pos in positions]
    ln = 0
    m = [0]
    n = 0
    for s in range(1, max(f)+1):
        t = [idx[0] + 1 for idx in enumerate(f) if idx[1] == s]
        l = len(t)
        for i in range(1, l+1):
            for j in range(t[i-1]-1, ln+i-1, -1):
                n = n+1
                if j == 1:
                    if n == 1:
                        zeros = [0]
                    else:
                        zeros = [0 for i in range(0, n)]
                    zeros[0:len(m)] = [m[i] for i in range(0, len(m))]
                    m = zeros
                    m[n-1] = f'lwij({k})'
                elif j == 2:
                    if n == 1:
                        zeros = [0]
                    else:
                        zeros = [0 for i in range(0, n)]
                    zeros[0:len(m)] = [m[i] for i in range(0, len(m))]
                    m = zeros
                    m[n-1] = f'(leye({logic})+lwij({k}))'
                else:
                    if n == 1:
                        zeros = [0]
                    else:
                        zeros = [0 for i in range(0, n)]
                    zeros[0:len(m)] = [m[i] for i in range(0, len(m))]
                    m = zeros
                    m[n-1] = f'(leye({k}^{j-1})+lwij({k}))'
                temp1 = f[j-1]
                f[j - 1] = f[j]
                f[j] = temp1
                temp2 = vars[j-1]
                vars[j-1] = vars[j]
                vars[j] = temp2
        ln += l
    k = len(f1)
    vars1[k:k+n] = m
    vars1[k + n: k + n + ln] = vars
    f1 = [item for sublist in [f1, [0 for i in range(0, n)], f] for item in sublist]
    return vars1, f1

def stdform(expr, options):
    """
    DESCRIPTION:
    A function transcribed from the MATLAB STP toolbox by Cheng Daizhan. At first we only assume the situation of a
    2-based logic.
    The variable logic sets the k-based logic, 2 for binary logic.
    """
    vars = expr.split(' ')
    f = [0 if var not in options else options.index(var) + 1 for var in vars]
    vars = ['MN', 'MD', 'C', 'F', 'A', 'B', 'MC', 'MC', 'MN', 'I', 'MN', 'C', 'MN', 'F', 'D', 'E', 'MN', 'MD', 'F', 'I',
            'G', 'H']
    f = [0, 0, 3, 6, 1, 2, 0, 0, 0, 9, 0, 3, 0, 6, 4, 5, 0, 0, 6, 9, 7, 8]
    # Set logic to 2 because we are in 2-based logic.
    vars, f = PR_Swap(vars, f, logic=2)
    vars, f = moveVarBack(vars, f, logic=2)
    vars, f = swapVarsByOrder(vars, f, logic=2)
    vars, f = PR_Swap(vars, f, logic=2)
    vars, f = moveVarBack(vars, f, logic=2)
    lm = convertSP(vars[0:len(f) - max(f)])
    positions = [u for u in range(0, len(f)) if u not in [
        c - 1 for c in range(0, 1 + len(f) - len([idx[0] + 1 for idx in enumerate(f) if idx[1] > 0]))
    ]]
    vars = [vars[pos] for pos in positions]
    return lm, vars

def stpn(matrices):
    """
    DESCRIPTION:
    A function to perform a left stp of n given matrices at the same time.
    :param matrices: [list] a list of matrices to be multiplied.
    :return: [numpy array] the result of the operations.
    """
    matrices.reverse()
    load = [matrices[0]]
    [load.append(stp(matrix, load[-1])) for matrix in matrices[1::]]
    return load[-1]

def matrix_eval(expr, variables):
    """
    DESCRIPTION:
    This function is devised to work in 2-based logic. That is the reason why we only gather a dictionary of pre-built
    functions.
    Dictionary:
    leye: logic identity matrix.
    MN: negation matrix.
    MC: conjunction matrix.
    MD: disjunction matrix.
    """
    # Identity matrix
    leye = np.array([[1, 0], [0, 1]])
    # Negation matrix
    MN = np.array([[0, 1], [1, 0]])
    # Conjunction matrix
    MC = np.array([[1, 0, 0, 0], [0, 1, 1, 1]])
    # Disjunction matrix
    MD = np.array([[1, 1, 1, 0], [0, 0, 0, 1]])

    return 3


def lGen(net, graph):
    """
    DESCRIPTION:
    A function to generate the L matrix from a given graph.
    :param net: [list] our representation of a network according with the formulation given in Murrugarra 2013.
    :param graph: [pandas DataFrame] the whole representation of the graph.
    :return: [numpy array] the L matrix of the graph.

    Dictionary:
    leye: logic identity matrix.
    MN: negation matrix.
    """
    # Iterate through the net
    for i in range(0, len(net)):
        # Set the structure of an INPUT node
        if net[i] == 'INPUT':
            net[i] = f'Meye {graph.index[i]}'
            continue
        # Iterate through the layers of a node
        load = ''
        steps = range(len(net[i])-1, -1, -1)
        m_step = min(steps)
        first_letter = list(net[i][m_step])
        for j in steps:
            # Create the structure string
            step = [item for sublist in [['MN', letter] for letter in list(net[i][j])] for item in sublist]
            ms = ['MN']
            if j == m_step and first_letter[0] not in graph.iloc[i]['activators']:
                ms = ['leye']
            step = ms + step
            load = ' '.join(step) + ' ' + load
        net[i] = load[0:len(load)-1]
    # Create the matrix from the string
    expr = ' '.join(net)
    options = list(graph.index)
    expr, variables = stdform(expr, options)
    l_matrix = matrix_eval(expr, variables)
    return l_matrix

