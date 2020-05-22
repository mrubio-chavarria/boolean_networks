

"""
INDICATIONS:
In this file, they are given several utilities to construct boolean networks based on the Semi-Tensor product. It is all
based in the structure stated in Murrugarra 2013.
"""
import re
import numpy as np
from numpy import matlib as mb
import numpy.matlib
from numpy import kron


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

def swapM(*args):
    """
    OLD BELOW IS THE BETTER ONE
    DESCRIPTION:
    A function that creates a swap matrix provided the dimensions for the indexes. Adapted from the one found in
    Analysis and Control of Boolean Networks.
    :param m: [int] first dimension.
    :param n: [int] second dimension.
    :return: [numpy array] the swap matrix.
    """

    # Parameters
    m = args[0]
    n = args[1] if len(args) > 1 else args[0]
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
                vars[v[j-2]-1] = convertKP(convertSP(vars[v[j-2]:v[j-1]-1]), logic)
                f[v[j-2]-1: v[j-2] + 1] = [0 for item in f[v[j-2]-1: v[j-2] + 1]]
                vars[v[j - 2]] = 'MR'
                k = range(v[j - 2] + 2, v[j-1], 1)  # BEWARE here
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
    load = matrices
    if len(matrices) > 1:
        matrices.reverse()
        load = [matrices[0]]
        [load.append(stp(matrix, load[-1])) for matrix in matrices[1::]]
    return load[-1]


def leye(n):
    """
    DESCRIPTION:
    A function to generate an identity matrix of nxn size.
    :param n: [int] dimensions number.
    :return: [numpy array] identity matrix.
    """
    return np.eye(n)

def leye(n):
    """
    DESCRIPTION:
    A function to generate an equivalence matrix of nxn size.
    :param n: [int] dimensions number.
    :return: [numpy array] identity matrix.
    """
    return np.eye(n)

def lwij(*args):
    """
    DESCRIPTION:
    A function that creates a swap matrix provided the dimensions for the indexes. Adapted from the one found in
    STP toolbox by Daizhan Cheng.
    :param args: [int] dimensions.
    :return: [numpy array] the swap matrix.
    """
    if len(args) == 1:
        m = args[0]
        n = m
    else:
        m = args[0]
        n = args[1]
    columns = np.eye(m*n)
    arrays = [np.array(range(1, m+1)) for k in range(0, n)]
    value = np.stack(arrays)
    value = value.flatten('F')
    i = np.transpose(value)
    j = numpy.matlib.repmat(np.array(range(1, n+1)), 1, m)[0, :]
    result = i + (j - 1)*m - 1
    arrays = [columns[:, i] for i in result]
    result = np.stack(arrays, 1)
    return result

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


def matrix_calc(expr, first_signal, second_signal, op_values, chunk_values, sym):
    r_pos = [i for i in range(0, len(expr)) if expr[i] == '-']
    s_pos = [i for i in range(0, len(expr)) if expr[i] == '+']
    sr_list = re.split('\+|\-', expr)
    results = []
    count = 0
    op_count = 0
    for cluster in sr_list:
        cluster = cluster.split('*')
        matrices = []
        for m in cluster:
            if m != second_signal and m != first_signal:
                matrices.append(sym[m])
            elif m == second_signal:
                matrices.append(chunk_values[count])
                count += 1
            else:
                matrices.append(op_values[op_count])
                op_count += 1

        results.append(stpn(matrices))
    result = results[0]
    for i in range(1, len(results)):
        m_r = min(r_pos) if len(r_pos) > 0 else max(s_pos) + 1
        m_s = min(s_pos) if len(s_pos) > 0 else max(r_pos) + 1
        if 0 < m_s < m_r:
            result = kron(result, results[i])
            s_pos.pop(s_pos.index(m_s))
        elif m_s > m_r > 0:
            raise ValueError
    return result


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

    # Identity matrix
    leye = np.array([[1, 0], [0, 1]])
    # Negation matrix
    MN = np.array([[0, 1], [1, 0]])
    # Conjunction matrix
    MC = np.array([[1, 0, 0, 0], [0, 1, 1, 1]])
    # Disjunction matrix
    MD = np.array([[1, 1, 1, 0], [0, 0, 0, 1]])
    # Equivalence matrix
    # Implication matrix
    # Power reducing matrix
    MR = np.array([[1, 0], [0, 0], [0, 0], [0, 1]])
    # XOR matrix

    # Dictionary
    sym = {'MN': MN, 'MC': MC, 'MD': MD, 'leye': leye, 'MR': MR}
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
    expr = ' '.join(net)
    options = list(graph.index)
    expr, variables = stdform(expr, options)
    l_matrix = matrix_eval(expr, variables)
    return l_matrix

