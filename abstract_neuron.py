from neuron import h
import os
import plotly
import matplotlib
import plotly.graph_objects as go
from abc import ABC, abstractmethod
import time
import numpy as np


class AbstractNeuron(ABC):
    """
    The abstract neuron class encapsulates all steps to initialize
    the model neuron. This includes loading the necessary files,
    inserting mechanisms, and defining init() for resetting the neuron.

    See the extensions of this class for models with different synapse
    types. Use only those models to run the simulations.

    NOTE: only one neuron per python session can be created.
    """

    def __init__(self):
        self.__hoc_dir = os.getcwd()
        self.create_cell()

    def create_cell(self):
        print("Creating cell")

        # Load standard run tools
        h.load_file("stdrun.hoc")

        # Load cell anatomical and biophysical properties
        self.load_hoc("nrnhoc/cellspec_c62564.hoc")
        self.load_hoc("nrnhoc/biophys.hoc")

        self.insert_mechanism("extracellular")
        self.insert_mechanism("xtrau")

        # Only interpolates sections that have xtrau
        self.load_hoc("nrnhoc/interpxyzu.hoc")

        # Automatically calls grindaway() in interpxyzu.hoc
        self.load_hoc("nrnhoc/setpointersu.hoc")

        # Computes scale factor used to calculate extracellular potential
        # produced by a uniform electrical field
        self.load_hoc("nrnhoc/calcrxcu.hoc")

        # Computes scale factor used to calculate extracellular potential
        # produced by a uniform electrical field
        self.load_hoc("nrnhoc/calcd.hoc")

        # Extracellular stimulus
        self.load_hoc("nrnhoc/zapstimu2.hoc")

        h('proc init() { nrnpython("AbstractNeuron._AbstractNeuron__neuron_init()") }')

        self.initialize()

    # noinspection PyMethodMayBeStatic
    def insert_mechanism(self, mechanism):
        for sec in h.allsec():
            sec.insert(mechanism)

    def load_hoc(self, hoc_file):
        hoc_path = os.path.join(self.__hoc_dir, hoc_file)
        h.load_file(hoc_path)

    @staticmethod
    def __neuron_init():
        h.t = 0

        for section in h.allsec():
            # Reset the resting potential
            section.v = h.Vrest

            # Set reversal potential for sodium channels
            if h.ismembrane("nax", sec = section) or \
               h.ismembrane("na3", sec = section):
                for segment in section: segment.ena = 55

            # Set reversal potential for potassium channels
            if h.ismembrane("kdr", sec = section) or \
               h.ismembrane("kap", sec = section) or \
               h.ismembrane("kad", sec = section):
                for segment in section: segment.ek = -90

            # Set reversal potential for h-current
            if h.ismembrane("hd", sec = section):
                for segment in section: segment.ehd_hd = -30

        # Set membrane potential to resting values
        h.finitialize(h.Vrest)
        # Calculate the currents
        h.fcurrent()

        for section in h.allsec():
            if h.ismembrane("na3", sec = section) or \
               h.ismembrane("nax", sec = section):
                for segment in section:
                    # Calculate passive current for sodium channels
                    segment.e_pas = segment.v + (segment.ina + segment.ik) / segment.g_pas

            if h.ismembrane("hd", sec = section):
                # Calculate passive current for h-current mechanisms
                for segment in section:
                    segment.e_pas = segment.e_pas + segment.i_hd / segment.g_pas

        print("Neuron initialized")

    """
        initialize() is called during create_cell(). Use this method to include
        synapses, set up stimulation parameters, etc.
    """
    @abstractmethod
    def initialize(self):
        pass

    def run(
            self,
            duration = 100,  # ms
            dt = 0.025  # ms
    ):
        self.__neuron_init()
        start = time.perf_counter()

        h.dt = dt
        h.tstop = duration
        h.continuerun(h.tstop)

        end = time.perf_counter()
        elapsed_time = end - start
        print(f"Elapsed time: {elapsed_time:.2f} seconds")

    @abstractmethod
    def get_synapse_info(self):
        pass

    def plot(self, draw_ef = False):
        # Do not show Burst cells or electric field components
        neuron_sections = h.SectionList([
            sec for sec in h.allsec() if "Burst" not in str(sec) and "sField" not in str(sec) and "sElec" not in str(sec)
        ])

        ps = h.PlotShape(neuron_sections, False)
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

        if draw_ef:
            #
            # Show direction of the electric field
            #
            x = h.sField.x3d(1)
            y = h.sField.y3d(1)
            z = h.sField.z3d(1)

            end = [x, y, z]
            start = [0, 0, 0]

            fig.add_trace(
                go.Cone(
                    x = [end[0]], y = [end[1]], z = [end[2]],
                    u = [end[0] - start[0]],
                    v = [end[1] - start[1]],
                    w = [end[2] - start[2]],
                    sizemode = "absolute",
                    sizeref = 50,
                    anchor = "tip",
                    colorscale = [[0, "red"], [1, "red"]],
                    showscale = False
                )
            )

            fig.add_trace(
                go.Scatter3d(
                    x = [start[0], end[0]],
                    y = [start[1], end[1]],
                    z = [start[2], end[2]],
                    mode = "lines",
                    line = dict(color = "red", width = 3)
                )
            )

            normal = [h.sField.x3d(1), h.sField.y3d(1), h.sField.z3d(1)]
            a, b, c = normal

            x = np.linspace(-1000, 1000, 10)
            y = np.linspace(-1000, 1000, 10)
            x, y = np.meshgrid(x, y)
            z = (-a * x - b * y) / c

            fig.add_trace(
                go.Surface(x = x, y = y, z = z, opacity = 0.25, colorscale = "gray", showscale = False)
            )

            # Clipping due to large plane size
            fig.update_layout(
                scene = dict(
                    zaxis = dict(
                        range = [-350, 350]
                    ),
                    yaxis = dict(
                        range = [-450, 550]
                    ),
                    xaxis = dict(
                        range = [-500, 500]
                    )
                )
            )

        synapse_info = self.get_synapse_info()
        if synapse_info is None or len(synapse_info) < 0:
            fig.show(config = { "scrollZoom": False })
            return

        for i, synapse in synapse_info.iterrows():
            dendrite_idx = int(synapse.Dendrite)
            dendrite_loc = synapse.Location

            i = int(h.apical_dendrite[dendrite_idx].n3d() * dendrite_loc)
            x = h.apical_dendrite[dendrite_idx].x3d(i)
            y = h.apical_dendrite[dendrite_idx].y3d(i)
            z = h.apical_dendrite[dendrite_idx].z3d(i)

            fig.add_trace(
                go.Scatter3d(
                    x = [x],
                    y = [y],
                    z = [z],
                    marker = dict(
                        color = "red",
                        size = 3
                    )
                )
            )

        fig.show(config = { "scrollZoom": False })
