

class InputAlterationException(Exception):
    """
    DESCRIPTION:
    The input to be launched when among the modified pathways, at least one of them has an input in its consequent.
    """
    # Methods
    def __init__(self):
        """
        DESCRIPTION:
        Builder of the class
        """
        super().__init__('Inputs cannot be altered')


class NotValidNetworkException(Exception):
    """
    DESCRIPTION:
    This exception is launched when in PyBoolNet it is not reached a coherent network given the expressions.
    """
    # Methods
    def __init__(self):
        """
        DESCRIPTION:
        Builder of the class
        """
        super().__init__('Not valid network found')


class ConvergenceError(Exception):
    """
    DESCRIPTION:
    This exception is launched when in the inference algorithm, or the conflicts resolution algorithm, does not reach
    a solution. It iterates beyond a limit.
    """
    # Methods
    def __init__(self):
        """
        DESCRIPTION:
        Builder of the class
        """
        super().__init__('Not solution found')


