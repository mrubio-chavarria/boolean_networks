#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
INDICATIONS:
In this script are encoded in Python all the classes described in MATLAB in Daizhan Cheng's STP toolbox. Available in
http://lsc.amss.ac.cn/~dcheng/
"""

import numpy as np
import math


class LM:
    """
    DESCRIPTION:
    The logical matrix class
    """
    def __init__(self, *args):
        """
        DESCRIPTION:
        Builder of the class. They can be introduced 1 or 2 arguments, if they are introduced more they will be ignored.
        They can be introduced a numpy logical matrix, or a numpy vector and an integer with the columns of the identity
        matrix.
        :param matrix: [numpy array] matrix filled with 1 or 0 to be converted in LM object.
        :parm v: [numpy array] vector with the columns in the identity matrix.
        :param n: [integer] number indicating the dimensions of the identity matrix.
        """
        # Check for LM objects
        ni = len(args)
        if ni and isinstance(args[0], LM):
            raise AttributeError('Use SET to modify the properties of LM objects.')
        self.v = []
        self.n = 0
        # Check for matrices
        if ni == 1:
            if ((args[0] == 1) | (args[0] == 0)).all() and\
                    ((np.sum(args[0], 0) == 1) | (np.sum(args[0], 0) == 0)).all():
                [p, q] = np.shape(args[0])
                k = np.where(np.transpose(args[0]).flatten() == 1)
                k = (np.array(np.where(np.transpose(args[0]).flatten() == 1)) + 1) % p
                self.v = (k + (k == 0)*p - 1)
                self.n = p
            else:
                raise AttributeError('They are valid logic matrices only.')
        elif ni > 1:
            if not isinstance(args[0], np.ndarray):
                raise ValueError('They are only allowed numpy arrays')
            if not isinstance(args[1], int):
                raise ValueError('They are only allowed integers')
            self.v = args[0]
            self.n = args[1]

    def __str__(self):
        """
        DESCRIPTION:
        Method to return a readable representation of the object.
        :return: [string] readable representation of the object.
        """
        return f'v: {self.v} n: {self.n}'

    def __mul__(self, other):
        """
        DECRIPTION:
        A method to overload the * operator according to the MATLAB code.
        :return: [LM] resul of the product.
        """
        # Prepare the case in which the other term is a logical matrix
        if not isinstance(other, LM):
            other = LM(other)
        # Parameters
        [m, n] = self.get_shape()
        [p, q] = other.get_shape()
        result = LM()
        # Check for incorrect matrices
        if n == 0 or q == 0:
            raise ValueError('Empty LM object')
        if n == p:
            result.v = self.v[other.v]
            result.n = self.n
        elif n % p == 0:
            k = int(n/p)
            r = np.transpose(np.reshape(self.v, [k, p]))
            result.v = np.transpose(r[:, other.v]).reshape([1, r[:, other.v].size])
            result.n = self.n
        elif p % n == 0:
            if n == 1:
                result.v = (self.v[0])*p + other.v
                result.n = p*m
            else:
                k = p / n
                x = np.array([math.floor(it) for it in (other.v/k).flatten()])
                t = (other.v + 1) % k
                y = t + (t == 0) * k
                result.v = self.v[0, x] * k + y - 1
                result.n = int(m * p / n)
        else:
            raise ValueError('Dimnensions must match multiple dimension condition')
        return result

    def __eq__(self, other):
        """
        DECRIPTION:
        A method to overload the == operator according to the MATLAB code.
        :return: [boolean] result of the comparison.
        """
        if isinstance(other, LM):
            if len(self.v) != len(other.v):
                raise ValueError('The column numbers must agree')
            result = self.v == other.v
        else:
            # To be IMPROVED
            raise ValueError('The two instances are to be of LM class')
        return result

    def get_shape(self):
        """
        DESCRIPTION:
        A method similar to the numpy shape, adapted to this object.
        :return: [tuple] dimensions of the matrix.
        """
        return self.n, len(self.v.flatten())



class STP:
    """
    DESCRIPTION:
    A class to hold a stp product in object form.
    """

