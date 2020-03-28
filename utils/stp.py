

"""
INDICATIONS:
In this file, they are given several utilities to construct boolean networks based on the Semi-Tensor product. It is all
based in the structure stated in Murrugarra 2013.
"""

import numpy as np
from numpy import matlib as mb


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

