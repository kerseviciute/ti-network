from neuron import h
import numpy as np
import time
from python.ca1 import NeuronCA1, Simulation

class Synapse:
  def __init__(
      self,
      neuron: NeuronCA1,
      dendrite: int,
      position: float,
      stimulus = None,
      delay: float = 0.0,
      initial_weight: float = 0.0,
      conductance: float = 0.0
    ):
    """
    Initialize a synapse.
    
    :param dendrite: dendrite index

    :param position: position along the dendrite, where the synapse is placed
    
    :param stimulus: Description
    
    :param delay: Description
    :type delay: float

    :param initial_weight: Description
    :type initial_weight: float

    :param conductance: Description
    :type conductance: float
    """
    self.neuron = neuron

    self.dendrite = dendrite
    self.position = position
    self.stimulus = stimulus
    self.delay = delay
    self.initial_weight = initial_weight
    self.conductance = conductance

    self.__initialize_synapses()

    # Setup to record what happens at the synapse
    self.__weights_0 = h.Vector().record(self.connection_ampa._ref_weight[0])
    self.__weights_1 = h.Vector().record(self.connection_ampa._ref_weight[1])
    self.__i_ampa = h.Vector().record(self.ampa._ref_i)
    self.__g_ampa = h.Vector().record(self.ampa._ref_g)
    self.__i_nmda = h.Vector().record(self.nmda._ref_i)
    self.__g_nmda = h.Vector().record(self.nmda._ref_g)
    self.__v = h.Vector().record(self.neuron.cell.apical_dendrite[dendrite](position)._ref_v)
    self.__t = h.Vector().record(h._ref_t)

    self.__stim_times = h.Vector()
    self.connection_ampa.record(self.__stim_times)

  def __initialize_synapses(self):
    self.ampa = h.STDPE2bis(self.position, sec = self.neuron.cell.apical_dendrite[self.dendrite])
    self.ampa.initial_weight = self.initial_weight * self.conductance
    self.nmda = h.nmdanet(self.position, sec = self.neuron.cell.apical_dendrite[self.dendrite])

    self.connection_ampa = self.__connect_to_stimulus(self.ampa)
    self.connection_nmda = self.__connect_to_stimulus(self.nmda)

  def __connect_to_stimulus(self, synapse):
    connection = h.NetCon(self.stimulus, synapse)
    connection.delay = self.delay
    connection.weight[0] = self.conductance

    return connection

  def set_conductance(self, conductance: float = 0.0):
    self.conductance = conductance
    self.connection_nmda.weight[0] = self.conductance
    self.connection_ampa.weight[0] = self.conductance

    self.ampa.initial_weight = self.initial_weight * self.conductance

  @property
  def weights(self):
    # weights_0 correspond to the maximum conductance
    # weights_1 correspond to the actual conductance
    # This returns the weight scaled from 0 to 1
    return np.array(self.__weights_1 / self.__weights_0)

  @property
  def time(self):
    return np.array(self.__t)

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

  @staticmethod
  def single_stimulus(time: float = 20.0):
    stimulus = h.NetStim()
    stimulus.number = 1
    stimulus.start = time
    stimulus.noise = 0
    stimulus.interval = 1e9

    return stimulus

  @staticmethod
  def generate_poisson_spike_times(seed, rate, duration):
    """Generate Poisson-distributed spike times."""
    np.random.seed(seed)

    spikes = []
    t = np.random.exponential(1000 / rate)
    spikes.append(t)

    while t < duration:
        spikes.append(t)
        isi = np.random.exponential(1000 / rate)  # Exponential ISI
        t += isi

    return spikes

  @staticmethod
  def poisson_stimulus(seed, rate, duration):
    spike_times = Synapse.generate_poisson_spike_times(seed, rate, duration)

    spike_vec = h.Vector(spike_times)
    stimulus = h.VecStim()
    stimulus.play(spike_vec)

    return stimulus

  def find_optimal_conductance(self, initial_conductance, target_epsp_diff = 0.15, epsilon = 0.025, verbose = False, run_verbose = False):
    if verbose:
      print(
        f"Estimating optimal conductance at the target EPSP {target_epsp_diff} mV with variability within {epsilon} " +
        f"mV (actual values between {round(target_epsp_diff - epsilon, 2)} and {round(target_epsp_diff + epsilon, 2)} mV)"
      )
    
    start = time.perf_counter()

    simulation = Simulation(neurons = [ self.neuron ])
    
    conductance = initial_conductance * 2
    epsp_diff = 0
    
    while abs(epsp_diff - target_epsp_diff) > epsilon:
      conductance_step = conductance / 2
      
      if epsp_diff > target_epsp_diff:
        conductance -= conductance_step
      else:
        conductance += conductance_step
      
      if verbose: print(f"Testing conductance {conductance} {conductance_step} {epsp_diff}")

      self.set_conductance(conductance = conductance)
      simulation.run(duration = 50, verbose = run_verbose)
      epsp_diff = np.abs(self.neuron.voltage[0] - self.neuron.voltage.max())
    
    end = time.perf_counter()
    elapsed_time = end - start
    if verbose:
      print(f"Elapsed time: {elapsed_time:.2f} seconds")
      print(f"Optimal conductance: {conductance} at EPSP {round(epsp_diff, 2)}")

    return { "conductance": conductance, "epsp": epsp_diff }
