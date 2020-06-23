#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
INDICATIONS:
In this script are encoded in Python all the algorithms described in MATLAB in Daizhan Cheng's STP toolbox. Available in
http://lsc.amss.ac.cn/~dcheng/
"""

import re
from math import prod

from pytictoc import TicToc

from stp.utils import lmc, lmn, lmr, lmd
from stp.classes import kron
from stp.utils import stpn, lwij, leye


def convertSP(expressions):
    if len(expressions) == 1:
        return expressions[0]
    else:
        chunk = '*'.join(expressions)
        return chunk


def convertKP(expressions, logic):
    return f'(leye({logic})+{expressions})'


def swap_vars_by_order(expressions, f, logic):
    k = logic
    v = [idx[0] + 1 for idx in enumerate(f) if idx[1] == 0]
    expressions1 = [expressions[i-1] for i in v]
    f1 = [f[i-1] for i in v]
    positions = [u for u in range(0, len(f)) if u not in [c - 1 for c in v]]
    expressions = [expressions[pos] for pos in positions]
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
                temp2 = expressions[j-1]
                expressions[j-1] = expressions[j]
                expressions[j] = temp2
        ln += l
    k = len(f1)
    expressions1[k:k+n] = m
    expressions1[k + n: k + n + ln] = expressions
    f1 = [item for sublist in [f1, [0 for i in range(0, n)], f] for item in sublist]
    return expressions1, f1


def pr_swap(expressions, f, logic):
    for i in range(1, max(f) + 1):
        v = [idx[0] + 1 for idx in enumerate(f) if idx[1] == i]
        j = 2
        while j <= len(v):
            chunk = f[v[j - 2]: v[j - 1] - 1]
            if v[j - 1] - v[j - 2] == 1:
                expressions[v[j - 2] - 1] = 'MR'
                f[v[j - 2] - 1] = 0
                v = [idx[0] + 1 for idx in enumerate(f) if idx[1] == i]
                continue
            elif not any(chunk):
                expressions[v[j - 2] - 1] = convertKP(convertSP(expressions[v[j - 2]:v[j - 1] - 1]), logic)
                f[v[j - 2] - 1: v[j - 2] + 1] = [0 for item in f[v[j - 2] - 1: v[j - 2] + 1]]
                expressions[v[j - 2]] = 'MR'
                k = range(v[j - 2] + 2, v[j - 1], 1)  # BEWARE here
                positions = [u for u in range(0, len(f)) if u not in [c - 1 for c in k]]
                f = [f[pos] for pos in positions]
                expressions = [expressions[pos] for pos in positions]
                v = [idx[0] + 1 for idx in enumerate(f) if idx[1] == i]
                continue
            j = j + 1
    return expressions, f


def swapvm(expressions, logic):
    value1 = convertSP(expressions[1::])
    value2 = convertKP(value1, logic)
    return [value2, expressions[0]]


def move_var_back(expressions, f, logic):
    ln = len(f)
    v = [idx[0] + 1 for idx in enumerate(f) if idx[1] > 0]
    l = len(v)
    for i in range(l, 0, -1):
        if v[i-1] == ln:
            continue
        chunk = [idx[0] + 1 for idx in enumerate(f[v[i-1]:ln-(l-i)]) if idx[1] == 0]
        if any(chunk):
            expressions[v[i-1]-1: v[i-1] + 1] = swapvm(expressions[v[i-1]-1: ln - (l - i)], logic)
            f[v[i-1]] = f[v[i-1]-1]
            f[v[i-1]-1] = 0
            if v[i-1] + 2 - (ln - (l - i)) == 0:
                k = [v[i-1] + 2]
            else:
                k = range(v[i-1] + 2, ln - (l - i) + 1, 1)
            if len(k) != 0:
                positions = [u for u in range(0, len(f)) if u not in [c-1 for c in k]]
                f = [f[pos] for pos in positions]
                expressions = [expressions[pos] for pos in positions]
            ln = len(f)
    return expressions, f


def stdform(expressions, options):
    """
    DESCRIPTION:
    A function transcribed from the MATLAB STP toolbox by Cheng Daizhan. At first we only assume the situation of a
    2-based logic. The variable logic sets the k-based logic, 2 for binary logic.
    :param expressions: [list] expressions of the different variables in string format within a list.
    :param options: [list] variables names.
    """

    # Parameters
    expressions = ' '.join(expressions).split(' ')
    f = [0 if var not in options else options.index(var) + 1 for var in expressions]

    # Algorithm
    # Set logic to 2 because we are in 2-based logic.
    expressions, f = pr_swap(expressions, f, logic=2)
    expressions, f = move_var_back(expressions, f, logic=2)
    expressions, f = swap_vars_by_order(expressions, f, logic=2)
    expressions, f = pr_swap(expressions, f, logic=2)
    expressions, f = move_var_back(expressions, f, logic=2)
    lm = convertSP(expressions[0:len(f) - max(f)])
    positions = [u for u in range(0, len(f)) if u not in [
        c - 1 for c in range(0, 1 + len(f) - len([idx[0] + 1 for idx in enumerate(f) if idx[1] > 0]))
    ]]
    variables = [expressions[pos] for pos in positions]
    return lm, variables


def detect_load(expr):
    loads_positions = []
    pos_open = [letter[0] for letter in enumerate(expr) if letter[1] == '(']
    pos_closed = [letter[0] for letter in enumerate(expr) if letter[1] == ')']
    for pos_c in pos_closed:
        pos_o = []
        i = 0
        while i < len(pos_open):
            pos = [pos_o for pos_o in pos_open if pos_o < pos_c]
            if len(pos) == 0:
                break
            pos = pos[-1]
            if pos_open[i] == pos:
                pos_o = pos_open.pop(i)
                i += -1
            i += 1
        loads_positions.append((pos_o, pos_c))
    return loads_positions


def parenthesis_filter(loads_positions, expr, tags):
    operations = []
    for pos in loads_positions[:]:
        tag = expr[pos[0] - 4:pos[0]]
        if tag in tags:
            new_pos = (pos[0] - 4, pos[1])  # Needed to reutilization of code
            # Ensure that we are not inside any $LOAD$
            condition = [True if pos[0] > l_pos[0] and pos[1] < l_pos[1] else False for l_pos in loads_positions]
            if not any(condition):
                condition = [True if pos[0] < l_pos[0] and pos[1] > l_pos[1] else False for l_pos in loads_positions]
                if not any(condition):
                    operations.append(new_pos)
    return operations


def expr_calc(expr):
    # Prepare the arguments
    args = expr[5:-1]
    if '^' in args:
        args = '**'.join(args.split('^'))
    args = eval(args)
    # Prepare the fucntion
    fun_name = expr[0:4]
    fun = globals()[fun_name]
    # Return the value
    value = fun(args)
    return value


def key2op(expr, operations, signal):
    op_values = [expr_calc(expr[pos[0]: pos[1] + 1]) for pos in operations]
    if len(operations) > 0:
        new_expr = expr[0:operations[0][0]] + signal
        if len(operations) == 1:
            new_expr = new_expr + expr[operations[0][1] + 1::]
        for i in range(1, len(operations)):
            first = expr[operations[i - 1][1] + 1: operations[i][0]]
            if i == len(operations) - 1:
                new_expr = new_expr + first + signal + expr[operations[i][1] + 1::]
                break
            new_expr = new_expr + first + signal
        expr = new_expr
    return expr, op_values


def embedded_filter(loads_positions):
    """
    DESCRIPTION:
    The aim of this function is to return the values of parenthesis which are not embedded in others parenthesis.
    """
    emb_poss = []
    for pos in loads_positions:
        for i in range(0, len(loads_positions)):
            i_pos = loads_positions[i]
            if i_pos[0] < pos[0] and i_pos[1] > pos[1]:
                emb_poss.append(pos)
    positions = [pos for pos in loads_positions if pos not in emb_poss]
    return positions


def matrix_calc(expr, first_signal, second_signal, op_values, chunk_values, sym):
    """
    DESCRIPTION:
    This function is the true eval within matrix eval.
    :param expr: [string] the expression to be assessed.
    :param first_signal: [string] mark for the position of values corresponding to functions.
    :param second_signal: [string] mark for the position of values corresponding to parenthesis.
    :param op_values: [list] matrices, the result of the functions.
    :param chunk_values: [list] matrices, the result of the parenthesis.
    :param sym: [dictionary] reference to values of precalculated matrices.
    """
    terms = []
    for term in expr.split(r'+'):
        factors = []
        for factor in term.split('*'):
            if factor == first_signal:
                factors.append(op_values.pop(0))
            elif factor == second_signal:
                factors.append(chunk_values.pop(0))
            else:
                factors.append(sym[factor])
        if len(factors) == 1:
            terms.append(factors[0])
        else:
            terms.append(prod(factors[1::], start=factors[0]))
    return sum(terms[1::], start=terms[0]) if len(terms) > 1 else terms[0]


def key2load(expr, positions, first_signal, second_signal):
    chunks = [expr[pos[0] + 1:pos[1]] for pos in positions]
    if len(positions) > 0:
        new_expr = expr[0:positions[0][0]] + second_signal
        if len(positions) == 1:
            new_expr = new_expr + expr[positions[0][1] + 1::]
        for i in range(1, len(positions)):
            first = expr[positions[i - 1][1] + 1: positions[i][0]]
            if i == len(positions) - 1:
                new_expr = new_expr + first + second_signal + expr[positions[i][1] + 1::]
                break
            new_expr = new_expr + first + second_signal
        expr = new_expr
    return expr, chunks


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

    Operations:
    +: addition
    -: subtraction
    *: multiplication (STP)
    """

    # Negation matrix
    MN = lmn()
    # Conjunction matrix
    MC = lmc()
    # Power reducing matrix
    MR = lmr()
    # Disjunction matrix
    MD = lmd()

    # Dictionary
    sym = {'MN': MN, 'MC': MC, 'MR': MR, 'MD': MD}
    # Tags
    tags = ['leye', 'lwij']
    # Catch the parenthesis by the keyword LOAD. Process the string
    loads_positions = detect_load(expr)
    # Filtering of function parenthesis
    operations = parenthesis_filter(loads_positions, expr, tags)
    # Change the parenthesis by the keyword OP
    first_signal = '$OP$'
    expr, op_values = key2op(expr, operations, first_signal)
    # Recalculate the positions
    loads_positions = detect_load(expr)
    # Filtering of embedded parenthesis
    positions = embedded_filter(loads_positions)
    # Make new expression with $LOAD$
    second_signal = '$LOAD$'
    expr, chunks = key2load(expr, positions, first_signal, second_signal)
    # Obtain the value of the parenthesis
    chunk_values = []
    for chunk in chunks:
        chunk_values.append(matrix_eval(chunk, variables))
    # Perform the operation
    result = matrix_calc(expr, first_signal, second_signal, op_values, chunk_values, sym)
    return result


def l_gen(net, graph, expressions=True):
    """
    DESCRIPTION:
    A function to generate the L matrix from a given graph.
    :param net: [list] our representation of a network according with the formulation given in Murrugarra 2013 or with a
    set of fixed expressions.
    :param graph: [pandas DataFrame] the whole representation of the graph.
    :param expressions: [boolean] a flag to indicate if the matrix is to be built from a set of expressions or from the
    structure of a NCBF.
    :return: [numpy array] the L matrix of the graph.
    """
    # Auxiliary function
    def aux_fun(net):
        for expr in net:
            expr = ''.join(r.findall(expr)).split(',')
            if expr[1] != '0' and expr[1] != '1':
                terms = list(aux_funII(expr[1].split('|')))
                expr = ' '.join(['MD ' + terms[i] if i != len(terms) - 1 else terms[i] for i in range(0, len(terms))])
            elif expr[1] == '0':
                expr = f'MC {graph.index[0]} MN {graph.index[0]}'
            elif expr[1] == '0':
                expr = f'MD {graph.index[0]} MN {graph.index[0]}'
            yield expr

    def aux_funII(terms):
        for term in terms:
            if term in list(graph.index):
                yield term
                continue
            # Probably this can be done faster without regex
            term = term if term[0] not in ['(', ')'] else term[1:-1]
            chunk = re.subn('[^\&]\w?', '', term)[1]
            header = '' if chunk == 1 and '!' in term else 'MC '
            term = re.sub('[\&](?=.*?\&)', '$MC$', term)
            term = re.sub('[$\&]', ' ', term.replace('!', 'MN$'))
            yield header + term

    # Choose a way to build the string
    r = re.compile(r'[^ ].*?')
    if expressions:
        net = list(aux_fun(net))
    else:
        # Iterate through the net
        for i in range(0, len(net)):
            # Set the structure of an INPUT node
            if net[i] == 'INPUT':
                net[i] = f'{graph.index[i]}'
                continue
            # Iterate through the layers of a node
            load = ''
            steps = range(len(net[i])-1, -1, -1)
            m_step = min(steps)
            first_letter = list(net[i][m_step])[0]
            for j in steps:
                layer = list(reversed(net[i][j]))
                # Assess if we are in the outer node or in the deepest
                ms = 'MN'
                if j == m_step and first_letter in graph.iloc[i]['activators']:
                    header = ms
                    if len(steps) != 1:
                        header += ' ' + 'MC'
                elif j == m_step and first_letter in graph.iloc[i]['inhibitors']:
                    header = ''
                    if len(steps) != 1:
                        header += 'MC'
                else:
                    if j == max(steps):
                        header = ms
                    else:
                        header = ms + ' ' + 'MC'
                # Calculate the whole load
                if len(layer) == 1:
                    level = 'MN' + ' ' + layer[0]
                else:
                    level = 'MN' + ' ' + layer[0]
                    for letter in layer[1::]:
                        level = 'MC' + ' ' + 'MN' + ' ' + letter + ' ' + level
                # Add inner load
                level = level + ' ' + load
                # Add the header
                load = header + ' ' + level
                pass
            node = load[0:len(load)-1].split(' ')
            if '' in node:
                node.pop(node.index(''))
            net[i] = ' '.join(node)
    # Create the matrix from the string
    options = list(graph.index)
    expr, variables = stdform(net, options)
    l_matrix = matrix_eval(expr, variables)
    return l_matrix

