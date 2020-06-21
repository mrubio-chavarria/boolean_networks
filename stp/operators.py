#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
INDICATIONS:
In this script are encoded in Python extra operators needed to perform with the classes. These algorithms are described
in MATLAB code in Daizhan Cheng's STP toolbox. Available in http://lsc.amss.ac.cn/~dcheng/
"""

from .classes import LM
import numpy as np

def kron(a, b):
    """
    DESCRIPTION:
    Kronecker product for LM class.
    :param a: [LM, numpy matrix] first matrix to be multiplied.
    :param b: [LM, numpy matrix] second matrix to be multiplied.
    :return: result of the multiplication.
    """
    # Check whether a or b are LM or not
    if isinstance(a, np.ndarray):
        a = LM(a)
    elif not isinstance(a, LM):
        raise ValueError('For Kronecker product, matrices are to be numpy array or LM.')
    if isinstance(b, np.ndarray):
        b = LM(b)
    elif not isinstance(b, LM):
        raise ValueError('For Kronecker product, matrices are to be numpy array or LM.')
    result = LM()
    # Obtain dimensions
    [m, n] = a.get_shape()
    [p, q] = b.get_shape()
    # Perform algorithm
    r = a.v*p
    r = np.ones([q, 1])*r
    t = np.transpose(b.v + 1)
    g = t * np.ones([1, n])
    r = r + g
    result.v = np.transpose(r).reshape((1, r.size)) - 1
    result.n = m*p
    return result


