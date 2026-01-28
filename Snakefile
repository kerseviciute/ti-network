import numpy as np

project = "ti-network"
neurons = np.arange(0, 5)

include: "rules/poisson.smk"

rule all:
  input:
    rules.poisson_10_1000_5_130_0.input,
    rules.poisson_10_0_0_0_0.input,
    synapse_info = expand("output/{project}/synapses/synapse_info.csv", project = project)

rule synapse_info:
  output:
    synapse_info = "output/{project}/synapses/synapse_info.csv",
    voltages = "output/{project}/synapses/synapse_voltage.csv",
    times = "output/{project}/synapses/times.csv"
  params:
    n_synapses = 200,
    initial_conductance = 0.5,
    target_epsp_diff = 0.15,
    epsilon = 0.025 # Around 0.15 with some variability
  conda: "neuron"
  script: "workflow/python/generate_synapses.py"

rule subthreshold_ef:
  output:
    subthreshold_ef = "output/{project}/subthreshold_ef/{carrier}_{beat}_{phi}_{psi}.csv"
  params:
    ti_carrier = "{carrier}",
    ti_beat = "{beat}",
    ti_phase = 10,
    ti_delay = 0,
    ti_phi = "{phi}",
    ti_psi = "{psi}",
    min_oscillations = 10,
    initial_amplitude = 2000,
    epsilon = 0.5
  conda: "neuron"
  script: "workflow/python/subthreshold_ef.py"
