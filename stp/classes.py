#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
INDICATIONS:
In this script are encoded in Python all the classes described in MATLAB in Daizhan Cheng's STP toolbox. Available in
http://lsc.amss.ac.cn/~dcheng/
"""

import copy
import numpy as np
import math


class LM:
    """
    DESCRIPTION:
    The logical matrix class
    """
    # Methods
    def __init__(self, *args):
        """
        DESCRIPTION:
        Constructor of the class. They can be introduced 1 or 2 arguments, if they are introduced more they will be ignored.
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
        self.v = np.array([[]]).astype(int)
        self.n = 0
        # Check for matrices
        if ni == 1:
            if ((args[0] == 1) | (args[0] == 0)).all() and (np.sum(args[0], 0) == 1).all():
                [p, q] = np.shape(args[0])
                k = (np.array(np.where(np.transpose(args[0]).flatten() == 1)) + 1) % p
                self.v = (k + (k == 0)*p - 1).astype(int)
                self.n = int(p)
            else:
                raise AttributeError('They are valid logic matrices only.')
        elif ni > 1:
            if not isinstance(args[0], np.ndarray):
                raise ValueError('They are only allowed numpy arrays')
            if not isinstance(args[1], int):
                raise ValueError('They are only allowed integers')
            self.v = np.reshape(args[0].astype(int), [1, args[0].size])
            self.n = int(args[1])

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
        # Check
        if not isinstance(other, LM):
            other = LM(other)
        self.v = np.reshape(self.v.astype(int), [1, self.v.size])
        other.v = np.reshape(other.v.astype(int), [1, other.v.size])
        # Parameters
        [m, n] = self.get_shape()
        [p, q] = other.get_shape()
        result = LM()
        # Check for incorrect matrices
        if n == 0 or q == 0:
            raise ValueError('Empty LM object')
        if n == p:
            try:
                result.v = self.v[0, other.v]
            except IndexError:
                print()
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
            raise ValueError('Dimensions must match multiple dimension condition')
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

    def __pow__(self, power, modulo=None):
        """
        DESCRIPTION:
        The power method.
        :param power: [integer] the power to which the object is to be powered.
        :return: [LM] the power of the object.
        """
        # Parameters
        [m, n] = self.get_shape()
        # Checks
        if not (m % n == 0 or n % m == 0):
            raise ValueError('The first input argument must match multiple dimension condition')
        elif not isinstance(power, int):
            raise ValueError('They power is to be an integer.')
        elif power < 0:
            raise ValueError('They power is to be a positive integer.')
        # Algorithm
        if power == 0:
            result = LM(np.array([[1]]))
        elif power == 1:
            result = self
        else:
            result = self * self ** (power - 1)
        return result

    def __add__(self, other):
        """
        DESCRIPTION:
        A method to overload the sum with the Kronecker product.
        :param other: [LM] matrix to be added.
        :return: [LM] addition of both matrices.
        """
        return kron(self, other)

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
    # Methods
    def __init__(self, *args):
        """
        DESCRIPTION:
        Constructor of the class.
        :param a: [numpy array, LM] matrix to be converted into a stp object.
        """
        self.c = np.array([])
        if args:
            # Check
            if not isinstance(args[0], np.ndarray):
                if not isinstance(args[0], LM):
                    raise ValueError('The STP constructor only admits numpy arrays. Invalid matrix type.')
                else:
                    self.c = args[0].v.astype(int)
            else:
                self.c = args[0].astype(int)

    def __str__(self):
        """
        DESCRIPTION:
        Method to return a readable representation of the object.
        :return: [string] readable representation of the object.
        """
        return f'c: {self.c}'

    def __eq__(self, other):
        """
        DESCRIPTION:
        Method to overload the == operator.
        :param other: [STP] the other STP object to be compared.
        :return: [numpy array] the matrix with the element-wise comparison.
        """
        # Check
        if not isinstance(other, STP):
            raise ValueError('Comparison can only be among STP objects. Invalid type.')
        return self.c == other.c

    def __sub__(self, other):
        """
        DESCRIPTION:
        Method to overload the - operator.
        :param other: [STP] the STP object to substract.
        :return: [numpy array] result of the substraction.
        """
        # Check
        if not isinstance(other, STP):
            raise ValueError('Substraction can only be among STP objects. Invalid type.')
        # Algorithm
        c = STP()
        if not (self.get_shape() == other.get_shape()):
            raise AttributeError('Inconsistent matrices sizes in STP substraction.')
        else:
            c.c = self.c - other.c
        return c

    def __mul__(self, other):
        """
        DESCRIPTION:
        Method to overload the * operator.
        :param other: [STP] the other factor of the product.
        :return: [STP] object resulting from the multiplication.
        """
        # Check
        if not isinstance(other, STP):
            raise ValueError('All the terms are to be STP to complete successfully the multiplication')
        [m, n] = self.get_shape()
        [p, q] = other. get_shape()
        if n == p:
            c = STP(self.c*other.c)
        else:
            t = np.lcm(n, p)
            c = STP(kron(self.c, np.eye(int(t/n)))*kron(other.c, np.eye(int(t/p))))
        return c

    def __add__(self, other):
        """
        DESCRIPTION:
        A method to overload the + operator.
        :param other: [STP, numpy array] the other factor to perform the addition.
        """
        # Check
        if not isinstance(other, STP):
            if isinstance(other, np.ndarray):
                other = STP(other)
            else:
                raise ValueError('Invalid type for addition. Valid types STP, numpy array and LM.')
        # Algorithm
        c = STP()
        if self.get_shape() == other.get_shape():
            c.c = self.c + other.c
        else:
            raise ValueError('Matrices dimensions must agree.')
        return c

    def __pow__(self, power, modulo=None):
        """
        DESCRIPTION:
        Method to overload the ** operator.
        :param power: [int] numeric exponent of the power.
        """
        # Check
        if not isinstance(power, int) or power < 0:
            raise ValueError('They are only allowed positive scalar in the exponent')
        # Algorithm
        if power == 0:
            result = STP(np.array([[1]]))
        elif power == 1:
            result = self
        else:
            result = self * self ** (power - 1)
        return result

    def get_shape(self):
        """
        DESCRIPTION:
        A method similar to the numpy shape, adapted to this object.
        :return: [tuple] dimensions of the matrix.
        """
        return np.shape(self.c)


"""
NOTES:
From here on they are encoded several operations that are needed for the methods of the classes above.
"""


def find(a, i):
    """
    DESCRIPTION:
    This function has not been taken from the MATLAB toolbox. This is a personal adaption of MATLAB 's find method.
    :param a: [numpy array] the vector to be assessed.
    :param i: [integer] the number to be found in a.
    :return: [numpy array] the vector with the selected positions.
    NOS HEMOS QUEDADO HACIENDO EL FIND PARA EL INV
    """


def inv(a):
    """
    DESCRIPTION:
    Function to calculate the inverse a logic matrix.
    :param a: [LM] logic matrix whose inverse is to be calculated.
    :return: [LM] inverse lm of a.
    """
    # Check
    if not issquare(a):
        raise ValueError('The input argument is not square.')
    if rank(a) != a.n:
        raise ValueError('The input argument is not invertible.')
    # Parameters
    result = LM()
    result.n = a.n
    # Algorithm
    result.v = np.array([np.where(a.v == i)[1] for i in range(0, result.n)]).reshape([1, result.n])
    return result


def rank(a):
    """
    DESCRIPTION:
    A function to determine the rank of a LM object.
    :param a: [LM] matrix to be analysed.
    :return: [boolean] indication of square or not.
    """
    return len(np.unique(a.v))


def unique(a):
    """
    DESCRIPTION:
    A function to return the unique elements of a in the form of another LM.
    :param a: [LM] matrix to be analysed.
    :return: [LM] matrix which the unique elements.
    """
    result = copy.deepcopy(a)
    result.v = np.unique(a.v)
    return result


def issquare(a):
    """
    DESCRIPTION:
    A function to determine whether a LM is square or not.
    :param a: [LM] matrix to be analysed.
    :return: [boolean] indication of square or not.
    """
    [m, n] = a.get_shape()
    return m == n


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
    a.v = a.v.reshape([1, a.v.size])
    b.v = b.v.reshape([1, b.v.size])
    # Obtain dimensions
    [m, n] = a.get_shape()
    [p, q] = b.get_shape()
    # Perform algorithm
    r = a.v*p
    r = np.ones([q, 1])*r
    t = np.transpose(b.v + 1)
    try:
        g = t * np.ones([1, n])
    except:
        print()
    r = r + g
    result.v = np.transpose(r).reshape((1, r.size)) - 1
    result.n = m*p
    return result


