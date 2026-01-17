from neuron import h
from dataclasses import dataclass
import plotly
import matplotlib
import time
import numpy as np


class Singleton(type):
  _instances = {}

  def __call__(cls, *args, **kwargs):
    if cls not in cls._instances:
      instance = super().__call__(*args, **kwargs)
      cls._instances[cls] = instance
    else:
      instance = cls._instances[cls]
      if hasattr(cls, '__allow_reinitialization') and cls.__allow_reinitialization:
        instance.__init__(*args, **kwargs)
    return instance


class NeuronLoader(metaclass = Singleton):
  """
    NeuronLoader class defines all files necessary
    to be loaded into a NEURON session. Since the files
    should be loaded only once during a session, the
    class is a singleton.
  """

  def __init__(self):
    # Load basic functions
    h.load_file("stdrun.hoc")

    # Only interpolates sections that have xtrau
    h.load_file("nrnhoc/interpxyzu.hoc")
    # Automatically calls grindaway() in interpxyzu.hoc
    h.load_file("nrnhoc/setpointersu.hoc")
    # Computes scale factor used to calculate extracellular potential
    # produced by a uniform electrical field
    h.load_file("nrnhoc/calcrxcu.hoc")
    # Computes scale factor used to calculate extracellular potential
    # produced by a uniform electrical field
    h.load_file("nrnhoc/calcd.hoc")
    # Extracellular stimulus
    h.load_file("nrnhoc/zapstimu2.hoc")

    # Load cell definition (contains morphological and biophysical properties)
    h.load_file("nrnhoc/geoc62564.hoc")


@dataclass
class TI:
  duration: float = 0
  delay: float = 0
  carrier: float = 0
  beat: float = 0
  amplitude: float = 0
  phase: float = 0
  psi_deg: float = 0
  phi_deg: float = 0
  __loader: NeuronLoader = NeuronLoader()

  def update_stimulus(self):
    h.setstim(
      self.delay,
      self.duration,
      self.carrier,
      self.carrier + self.beat,
      self.amplitude,
      self.phase
    )

  def update_field(self):
    field_length = 100
    h.changefield(field_length, self.psi_deg, self.phi_deg)

  def __setattr__(self, name, value):
    super().__setattr__(name, value)
    if name in {
      "duration", "delay", "carrier", "beat", "amplitude", "phase", "psi_deg", "phi_deg"
    }:
      if getattr(self, "_ti_ready", False):
        self.update_stimulus()
        self.update_field()

  def __post_init__(self):
    object.__setattr__(self, "_ti_ready", True)
    self.update_stimulus()
    self.update_field()


class NeuronCA1:
  def __init__(self, ap_threshold = -20):
    self.loader = NeuronLoader()

    # Load the cell with morphological and biophysical properties
    self.cell = h.CA1()
    self.insert_mechanism("extracellular")
    self.insert_mechanism("xtrau")

    h.setpointers()

    self.__v = h.Vector().record(self.cell.soma[0](0.5)._ref_v)
    self.__t = h.Vector().record(h._ref_t)

    # Record spike times
    self.__spike_times = h.Vector()
    ap_connection = h.NetCon(self.cell.soma[0](0.5)._ref_v, None, sec = self.cell.soma[0])
    ap_connection.threshold = ap_threshold
    ap_connection.record(self.__spike_times)

  @property
  def voltage(self):
    return np.array(self.__v)

  @property
  def time(self):
    return np.array(self.__t)

  @property
  def spike_times(self):
    if len(self.__spike_times) > 0:
      return np.array(self.__spike_times)
    else:
      return np.array([])

  def insert_mechanism(self, mechanism):
    for section in self.cell.all:
      section.insert(mechanism)

  def reset_state(self):
    h.t = 0

    for section in self.cell.all:
      # Reset the resting potential
      section.v = self.cell.Vrest

      # Set reversal potential for sodium channels
      if h.ismembrane("nax", sec = section) or h.ismembrane("na3", sec = section):
        for segment in section: segment.ena = 55

      # Set reversal potential for potassium channels
      if h.ismembrane("kdr", sec = section) or h.ismembrane("kap", sec = section) or h.ismembrane("kad", sec = section):
        for segment in section: segment.ek = -90

      # Set reversal potential for h-current
      if h.ismembrane("hd", sec = section):
        for segment in section: segment.ehd_hd = -30

    # Set membrane potential to resting values
    h.finitialize(self.cell.Vrest)
    # Calculate the currents
    h.fcurrent()

    for section in self.cell.all:
      if h.ismembrane("na3", sec = section) or h.ismembrane("nax", sec = section):
        for segment in section:
          # Calculate passive current for sodium channels
          segment.e_pas = segment.v + (segment.ina + segment.ik) / segment.g_pas

      if h.ismembrane("hd", sec = section):
        # Calculate passive current for h-current mechanisms
        for segment in section:
          segment.e_pas = segment.e_pas + segment.i_hd / segment.g_pas

    print("Neuron state has been reset")

  def translate(self, dx = 0.0, dy = 0.0, dz = 0.0):
    """
    Translate all pt3d points of every section in `cell` by (dx, dy, dz).
    Works only if sections have 3D points (pt3dadd).
    """
    for sec in self.cell.all:
      n = int(h.n3d(sec = sec))
      if n == 0:
        continue

      for i in range(n):
        x = h.x3d(i, sec = sec) + dx
        y = h.y3d(i, sec = sec) + dy
        z = h.z3d(i, sec = sec) + dz
        d = h.diam3d(i, sec = sec)
        h.pt3dchange(i, x, y, z, d, sec = sec)

  def plot(self):
    sections = h.SectionList([sec for sec in self.cell.all])
    NeuronCA1.__plot(sections)

  @staticmethod
  def plot_all():
    sections = h.SectionList([sec for sec in h.allsec()])
    NeuronCA1.__plot(sections)

  @staticmethod
  def __plot(sections):
    ps = h.PlotShape(sections, False)
    ps.show(1)

    fig = ps.plot(plotly, cmap = matplotlib.colormaps["inferno"])

    fig.update_layout(
      scene = dict(
        xaxis_title = "",
        yaxis_title = "",
        zaxis_title = "",
        xaxis = dict(showbackground = False, showticklabels = False),
        yaxis = dict(showbackground = False, showticklabels = False),
        zaxis = dict(showbackground = False, showticklabels = False)
      )
    )

    fig.show()


class Simulation:
  def __init__(self, neurons, ti = TI()):
    self.neurons = neurons
    self.ti = ti

  def get_dt(self):
    return 0.025 if self.ti.carrier < 5000 else 0.0025

  def run(self, duration):
    for neuron in self.neurons:
      neuron.reset_state()

    start = time.perf_counter()

    h.dt = self.get_dt()
    h.tstop = duration
    h.run(h.tstop)

    end = time.perf_counter()
    elapsed_time = end - start
    print(f"Elapsed time: {elapsed_time:.2f} seconds")
