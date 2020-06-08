
class Result:
    """
    DESCRIPTION:
    A class built to store the results of the simulation.
    """

    def __init__(self, network, pathways, maps_set, conflicts, simulation, iterations):
        """
        DESCRIPTION:
        Constructor to the Result class.
        :param network: [list] the structured network employed.
        :param pathways: [list] objects pathway, the result of the simulation.
        :param conflicts: [list] objects conflict, the record of the simulation.
        :param maps_set: [KMap] the set of maps employed during the simulation.
        :param simulation: [int] identifier of the simulation.
        :parm iterations: [int] number of the iterations performed.
        """
        self.network = network
        self.pathways = None
        self.set_pathways(pathways)
        self.maps_set = maps_set
        self.conflicts = conflicts
        self.code = None
        self.code = self.get_code()
        self.iterations = iterations
        self.simulation = simulation

    def get_code(self):
        """
        DESCRIPTION:
        A method to built a unique identifier of the type to which the result belongs.
        :return: [string] code of the kind of pathway
        """
        if self.code is None:
            net = ''.join([it for sl in self.network for it in sl])
            codes = ['$']
            [codes.append(f'{str(hex(int("".join(pathway.region_of_interest), 2)))}:{pathway.consequent}|')
             for pathway in self.pathways]
            code = ''.join(codes) + f'{net}|$'
        else:
            code = self.code
        return code

    def set_pathways(self, pathways):
        """
        DESCRIPTION:
        A method to get the pathways always ordered in the same manner.
        :param pathways: [list] objects pathway, the result of the simulation.
        """
        pathways.sort(key=lambda x: x.consequent)
        self.pathways = pathways

