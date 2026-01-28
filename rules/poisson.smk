rule poisson:
  output:
    neuron_voltage = "output/{project}/poisson/{rate}_{carrier}_{beat}_{amplitude}_{weight}/{run_id}_neuron_voltage.csv",
    neuron_spikes = "output/{project}/poisson/{rate}_{carrier}_{beat}_{amplitude}_{weight}/{run_id}_neuron_spikes.csv",
    synapse_weights = "output/{project}/poisson/{rate}_{carrier}_{beat}_{amplitude}_{weight}/{run_id}_synapse_weights.csv",
    synapse_voltages = "output/{project}/poisson/{rate}_{carrier}_{beat}_{amplitude}_{weight}/{run_id}_synapse_voltages.csv",
    synapse_stimuli = "output/{project}/poisson/{rate}_{carrier}_{beat}_{amplitude}_{weight}/{run_id}_synapse_stimuli.csv"
  params:
    n_synapses = 10,
    run_id = "{run_id}",
    carrier = "{carrier}",
    beat = "{beat}",
    phase = 10,
    delay = 0,
    phi_deg = 90,
    psi_deg = 90,
    amplitude = "{amplitude}",
    duration = 1000,
    rate = "{rate}",
    weight = "{weight}"
  conda: "neuron"
  script: "../workflow/python/poisson.py"

rule poisson_10_1000_5_130_0:
  input:
    expand(
      "output/{project}/poisson/{rate}_{carrier}_{beat}_{amplitude}_{weight}/{run_id}_neuron_voltage.csv",
      project = project,
      rate = 10,
      carrier = 1000,
      beat = 5,
      amplitude = 130,
      weight = 0,
      run_id = neurons
    )

rule poisson_10_0_0_0_0:
  input:
    expand(
      "output/{project}/poisson/{rate}_{carrier}_{beat}_{amplitude}_{weight}/{run_id}_neuron_voltage.csv",
      project = project,
      rate = 10,
      carrier = 0,
      beat = 0,
      amplitude = 0,
      weight = 0,
      run_id = neurons
    )
