from neuron import h
import numpy as np


class SynapseCA3:
    # NOTE: the synaptic model used here cannot undergo LTD!

    def __init__(
            self,
            dendrite_loc,
            dendrite_idx,
            stimulus,
            delay: float = 0.0,
            initial_weight: float = 0.0,
            initial_conductance: float = 0.0
    ):
        self.dendrite = dendrite_idx
        self.position = dendrite_loc

        self.initial_weight = initial_weight

        self.ampa = h.STDPE2bis(dendrite_loc, sec = h.apical_dendrite[dendrite_idx])
        self.ampa.initial_weight = self.initial_weight * initial_conductance
        self.nmda = h.nmdanet(dendrite_loc, sec = h.apical_dendrite[dendrite_idx])
        self.stimulus = stimulus
        self.delay = delay

        self.connection_ampa = self.__connect(
            self.ampa, self.stimulus, delay = delay, initial_conductance = initial_conductance
        )
        # TODO: what should the threshold be?
        self.connection_nmda = self.__connect(
            self.nmda, self.stimulus, delay = delay, initial_conductance = initial_conductance
        )

        self.__weights_0 = h.Vector().record(self.connection_ampa._ref_weight[0])
        self.__weights_1 = h.Vector().record(self.connection_ampa._ref_weight[1])
        self.__i_ampa = h.Vector().record(self.ampa._ref_i)
        self.__g_ampa = h.Vector().record(self.ampa._ref_g)
        self.__i_nmda = h.Vector().record(self.nmda._ref_i)
        self.__g_nmda = h.Vector().record(self.nmda._ref_g)
        self.__v = h.Vector().record(h.apical_dendrite[dendrite_idx](dendrite_loc)._ref_v)

        self.__stim_times = h.Vector()
        self.connection_ampa.record(self.__stim_times)

    def set_initial_conductance(self, conductance: float):
        self.connection_nmda.weight[0] = conductance
        self.connection_ampa.weight[0] = conductance

        self.ampa.initial_weight = self.initial_weight * conductance

    @staticmethod
    def __connect(synapse, stimulus, delay: float = 0.0, initial_conductance: float = 1.0):
        connection = h.NetCon(stimulus, synapse)
        connection.delay = delay
        connection.weight[0] = initial_conductance

        return connection

    @property
    def weights(self):
        # weights_0 correspond to the maximum conductance
        # weights_1 correspond to the actual conductance
        return np.array(self.__weights_1 / self.__weights_0)

    @property
    def ampa_current(self):
        return np.array(self.__i_ampa)

    @property
    def ampa_conductance(self):
        return np.array(self.__g_ampa)

    @property
    def voltage(self):
        return np.array(self.__v)

    @property
    def nmda_current(self):
        return np.array(self.__i_nmda)

    @property
    def nmda_conductance(self):
        return np.array(self.__g_nmda)

    @property
    def stimuli(self):
        if len(self.__stim_times) > 0:
            return np.array(self.__stim_times) + self.delay
        else:
            return np.array([])
