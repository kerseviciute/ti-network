import matplotlib.pyplot as plt
import numpy as np


class NeuronPlots:
    @staticmethod
    def plot_all_synapses(neuron):
        fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (10, 10))
        axs = axs.flatten()

        for synapse in neuron.synapses:
            axs[0].plot(neuron.time, synapse.weights, linewidth = 0.5)

        axs[0].set_ylim(-0.1, 1.1)
        axs[0].set_title("Synaptic weights")

        for synapse in neuron.synapses:
            axs[1].plot(neuron.time, synapse.ampa_current, linewidth = 0.5)

        axs[1].set_title("Synapse current")

        for synapse in neuron.synapses:
            axs[2].plot(neuron.time, synapse.ampa_conductance, linewidth = 0.5)

        axs[2].set_title("Synapse conductance")

        for synapse in neuron.synapses:
            axs[3].plot(neuron.time, synapse.voltage, linewidth = 0.5)

        axs[3].set_title("Synapse EPSP")

    @staticmethod
    def plot_all_nmda(neuron):
        fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (10, 5))
        axs = axs.flatten()

        for synapse in neuron.synapses:
            axs[0].plot(neuron.time, synapse.nmda_current, linewidth = 0.5)

        axs[0].set_title("Synapse current")

        for synapse in neuron.synapses:
            axs[1].plot(neuron.time, synapse.nmda_conductance * 1000, linewidth = 0.5)

        axs[1].set_xlabel("Time, ms")
        axs[1].set_ylabel(r"Conductance, nS")
        axs[1].set_title("Synapse conductance")

        fig.tight_layout()

    @staticmethod
    def plot_conductance(neuron):
        fig, axs = plt.subplots(nrows = 1, ncols = 3, figsize = (12, 4))
        axs = axs.flatten()

        for synapse in neuron.synapses:
            axs[0].plot(neuron.time, synapse.ampa_conductance * 1000, linewidth = 0.5)

        axs[0].set_xlabel("Time, ms")
        axs[0].set_ylabel(r"Conductance, nS")
        axs[0].set_title("AMPA conductance")

        for synapse in neuron.synapses:
            axs[1].plot(neuron.time, synapse.nmda_conductance * 1000, linewidth = 0.5)

        axs[1].set_xlabel("Time, ms")
        axs[1].set_ylabel(r"Conductance, nS")
        axs[1].set_title("NMDA conductance")

        for synapse in neuron.synapses:
            axs[2].plot(neuron.time, (synapse.ampa_conductance + synapse.nmda_conductance) * 1000, linewidth = 0.5)

        axs[2].set_xlabel("Time, ms")
        axs[2].set_ylabel(r"Conductance, nS")
        axs[2].set_title("AMPA and NMDA conductance")

        fig.tight_layout()

    @staticmethod
    def plot_soma_voltage(neuron, scaled = True, ax = None):
        if ax is None:
            fig, ax = plt.subplots()

        ax.plot(neuron.time, neuron.voltage, linewidth = 0.5, color = "black")
        ax.set_xlabel("Time, ms")
        ax.set_ylabel("Voltage, mV")
        ax.set_xlim(0, neuron.time.max())

        if scaled:
            ax.set_ylim(-100, 70)

        return ax

    @staticmethod
    def plot_spike_frequency(neuron):
        inst_freq, time_points = neuron.spike_frequency

        fig, axs = plt.subplots(nrows = 2, ncols = 1, figsize = (8, 5))

        axs[0].plot(neuron.time, neuron.voltage, linewidth = 0.5, color = "black")
        axs[0].set_xlabel("Time (ms)")
        axs[0].set_ylabel("Voltage (mV)")
        axs[0].set_xlim(np.min(neuron.time), np.max(neuron.time))
        axs[0].set_ylim(-100, 70)

        axs[1].scatter(time_points, inst_freq, s = 3, color = "black")
        axs[1].set_xlabel("Time (ms)")
        axs[1].set_ylabel("Frequency (Hz)")
        axs[1].set_xlim(np.min(neuron.time), np.max(neuron.time))

        if len(inst_freq) > 0:
            axs[1].set_ylim(-1, np.max(inst_freq) + 5)

        fig.tight_layout()

    @staticmethod
    def plot_synaptic_weights(neuron):
        for synapse in neuron.synapses:
            plt.plot(neuron.time, synapse.weights)

        plt.ylim(-0.1, 1.1)
        plt.title("Synaptic weights")
        plt.xlabel("Time, ms")
        plt.ylabel("Synaptic weight")
