# Include path to neuron classes
import subprocess; git_root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text = True).strip()
import sys; sys.path.append(git_root)

from python.ca1 import NeuronCA1, Simulation
from python.synapse import Synapse
import pandas as pd
import random
import numpy as np

n_synapses = int(snakemake.params["n_synapses"])
run_id = int(snakemake.params["run_id"])
beat = int(snakemake.params["beat"])
carrier = int(snakemake.params["carrier"])
phase = int(snakemake.params["phase"])
delay = int(snakemake.params["delay"])
phi_deg = int(snakemake.params["phi_deg"])
psi_deg = int(snakemake.params["psi_deg"])
amplitude = int(snakemake.params["amplitude"])
duration = int(snakemake.params["duration"])
rate = int(snakemake.params["rate"])
weight = float(snakemake.params["weight"])

# Set random seed as this run's ID
np.random.seed(run_id)
stimulus_seeds = np.random.choice(range(10_000), n_synapses, replace = False)

synapse_info = pd.read_csv(snakemake.input["synapse_info"], index_col = 0)
synapse_info = synapse_info.iloc[ (n_synapses * run_id):(n_synapses * (run_id + 1)) ].reset_index(drop = True)

# Create neuron and insert all synapses
neuron = NeuronCA1()
synapses = [ 
  Synapse(
    neuron = neuron,
    dendrite = int(info.dendrite),
    position = info.position,
    conductance = info.conductance,
    initial_weight = weight,
    stimulus = Synapse.poisson_stimulus(seed = stimulus_seeds[i], rate = rate, duration = duration)
  ) for i, info in synapse_info.iterrows()
]

# Set up simulation
simulation = Simulation([ neuron ])

simulation.ti.beat = beat
simulation.ti.carrier = carrier
simulation.ti.phase = phase
simulation.ti.delay = delay
simulation.ti.duration = duration
simulation.ti.phi_deg = phi_deg
simulation.ti.psi_deg = psi_deg
simulation.ti.amplitude = amplitude

simulation.run(duration = duration)

################ Save the data

pd.DataFrame({
  "time": neuron.time,
  "voltage": neuron.voltage,
}).to_csv(snakemake.output["neuron_voltage"])

pd.DataFrame({
  "time": neuron.spike_times
}).to_csv(snakemake.output["neuron_spikes"])

# Collect all synapse weights
weights = []
for synapse in synapses:
  weights.append(synapse.weights)
weights = np.array(weights)

# Since weights are mostly stable, keep only unique values
weights, index = np.unique(weights, axis = 1, return_index = True)
time = synapses[0].time[index].reshape(-1)
weights = np.vstack([time, weights])

np.savetxt(
  snakemake.output["synapse_weights"],
  weights,
  delimiter = ",",
  fmt = "%g"
)

stimuli = []
for synapse in synapses:
    stimuli.append(synapse.stimuli)

stimuli = [list(row) for row in stimuli]
max_len = max(len(row) for row in stimuli)
padded_data = [row + [""] * (max_len - len(row)) for row in stimuli]

np.savetxt(
  snakemake.output["synapse_stimuli"],
  padded_data,
  delimiter = ",",
  fmt = "%s"
)

