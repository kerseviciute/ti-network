# Include path to neuron classes
import subprocess; git_root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text = True).strip()
import sys; sys.path.append(git_root)

from python.ca1 import NeuronCA1
from python.synapse import Synapse
import pandas as pd
import random
import numpy as np

n_synapses = int(snakemake.params["n_synapses"])
initial_conductance = float(snakemake.params["initial_conductance"])
target_epsp_diff = float(snakemake.params["target_epsp_diff"])
epsilon = float(snakemake.params["epsilon"])

neuron = NeuronCA1()

random.seed(23875934)
synapse_locations = pd.concat(
  [pd.DataFrame([neuron.generate_synapse_location()]) for _ in range(n_synapses)],
  ignore_index = True
)
synapse_locations["conductance"] = 0.0
synapse_locations["epsp"] = 0.0

voltages = []
times = []
for i, synapse_loc in synapse_locations.iterrows():
  print(synapse_loc)
  
  synapse = Synapse(
    neuron = neuron,
    dendrite = int(synapse_loc.dendrite),
    position = synapse_loc.position,
    stimulus = Synapse.single_stimulus()
  )

  estimate = synapse.find_optimal_conductance(
    initial_conductance = initial_conductance,
    target_epsp_diff = target_epsp_diff,
    epsilon = epsilon,
    verbose = True
  )

  synapse_locations.loc[synapse_locations.index[i], "conductance"] = estimate["conductance"]
  synapse_locations.loc[synapse_locations.index[i], "epsp"] = estimate["epsp"]

  voltages.append(synapse.voltage.copy())
  times.append(synapse.time.copy())

  # Make sure synapse is deleted from neuron (to not influence optimal 
  # conductance estimation for other synapses)
  del synapse

np.savetxt(
  snakemake.output["voltages"],
  np.array(voltages),
  delimiter = ",",
  fmt = "%g"
)

np.savetxt(
  snakemake.output["times"],
  np.array(times),
  delimiter = ",",
  fmt = "%g"
)

synapse_locations.to_csv(snakemake.output["synapse_info"])
