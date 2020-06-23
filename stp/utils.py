#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
INDICATIONS:
In this file are the functions needed to generate the different kinds of matrices employed in the algorithms. Adapted
from http://lsc.amss.ac.cn/~dcheng/
"""

from stp.classes import LM, kron
import numpy as np


def stp(a, b):
    """
    DESCRIPTION:
    A function to perform an stp product between two matrices.
    :param a: [numpy array] the first matrix to be multiplied.
    :param b: [numpy array] the second matrix to be multiplied.
    """
    # Check
    if not isinstance(a, LM):
        if not isinstance(a, np.ndarray):
            raise ValueError('Invalid data type for STP. Only np.array or LM are allowed.')
        a = LM(a)
    if not isinstance(b, LM):
        if not isinstance(b, np.ndarray):
            raise ValueError('Invalid data type for STP. Only np.array or LM are allowed.')
        b = LM(b)
    [m, n] = a.get_shape()
    [p, q] = b.get_shape()
    if n == p:
        c = a * b
    else:
        z = np.lcm(n, p)
        c = kron(a, leye(int(z/n)))*kron(b, leye(int(z/p)))
    return c


def stpn(matrices):
    """
    DESCRIPTION:
    A function to perform a left stp of n given matrices at the same time.
    :param matrices: [list] a list of matrices to be multiplied.
    :return: [numpy array] the result of the operations.
    """
    load = matrices
    if len(matrices) > 1:
        matrices.reverse()
        load = [matrices[0]]
        [load.append(stp(matrix, load[-1])) for matrix in matrices[1::]]
    return load[-1]


def lmd(k=None):
    """
    DESCRIPTION:
    A function to generate disjunction matrices of the class LM.
    :param k: [int] to set the k-valued logic.
    :return: [LM] disjunction matrix.
    """
    # Check
    k = 2 if k is None else k
    # Algorithm
    result = LM()
    result.n = k
    a = range(1, k+1)
    p = np.repeat(a, repeats=k, axis=0)
    q = np.tile(a, reps=k)
    b = p <= q
    result.v = np.multiply(p, b) + np.multiply(q, 1 - b) - 1
    result.v = np.reshape(result.v, [1, result.v.size])
    return result


def lmc(k=None):
    """
    DESCRIPTION:
    A function to generate conjunction matrices of the class LM.
    :param k: [int] to set the k-valued logic.
    :return: [LM] conjunction matrix.
    """
    # Check
    k = 2 if k is None else k
    # Algorithm
    m = LM()
    m.n = k
    a = range(1, k+1)
    p = np.repeat(a, repeats=k, axis=0)
    q = np.tile(a, reps=k)
    b = (p >= q).astype(int)
    m.v = np.multiply(p, b) + np.multiply(q, 1 - b) - 1
    m.v = np.reshape(m.v, [1, m.v.size])
    return m


def lmn(k=None):
    """
    DESCRIPTION:
    A function to generate negation matrices of the class LM.
    :param k: [int] to set the k-valued logic.
    :return: [LM] negation matrix.
    """
    # Check
    k = 2 if k is None else k
    return LM(np.array([range(k, 0, -1)]) - 1, k)


def lwij(*args):
    """
    DESCRIPTION:
    A function to generate swap matrices of the class LM.
    :param args[0] (M): [int] first dimension.
    :param args[1] (N): [int] second dimension.
    :return: [LM] swap matrix of MNxMN dimensions.
    """
    # Check
    if len(args) == 1:
        m = n = args[0]
    else:
        m = args[0]
        n = args[1]
    # Algorithm
    w = LM()
    w.n = m*n
    i = np.repeat(range(1, m+1), repeats=n, axis=0)
    j = np.tile(range(1, n+1), reps=m)
    w.v = i + (j-1)*m - 1
    w.v = np.reshape(w.v, [1, w.v.size])
    return w


def lmr(k=None):
    """
    DESCRIPTION:
    A function to generate power reducing matrices of the class LM.
    :param k: [int] to set the k-valued logic.
    :return: [LM] power reducing matrix.
    """
    # Check
    k = 2 if k is None else k
    # Algorithm
    a = np.array([range(1, k+1)])
    mr = LM(a+(a-1)*k-1, k**2)
    return mr


def leye(k=None):
    """
    DESCRIPTION:
    A function to generate identity matrices of the class LM.
    :param k: [int] to set the k-valued logic.
    :return: [LM] identity matrix.
    """
    # Check
    k = 2 if k is None else k
    return LM(np.array([range(1, k+1)]) - 1, k)


