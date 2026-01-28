# Include path to neuron classes
import subprocess; git_root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text = True).strip()
import sys; sys.path.append(git_root)

from python.ca1 import NeuronCA1, Simulation
import pandas as pd

neuron = NeuronCA1()
simulation = Simulation(neurons = [ neuron ])

simulation.ti.beat = float(snakemake.params["ti_beat"])
simulation.ti.carrier = int(snakemake.params["ti_carrier"])
simulation.ti.phase = int(snakemake.params["ti_phase"])
simulation.ti.delay = int(snakemake.params["ti_delay"])

simulation.ti.phi_deg = int(snakemake.params["ti_phi"])
simulation.ti.psi_deg = int(snakemake.params["ti_psi"])

min_oscillations = int(snakemake.params["min_oscillations"])

initial_amplitude = int(snakemake.params["initial_amplitude"])
epsilon = float(snakemake.params["epsilon"])

duration = simulation.ti.get_optimal_duration(min_oscillations = min_oscillations)
simulation.ti.duration = duration
print(f"Estimated optimal simulation duration: {duration} ms")

amplitude = initial_amplitude
amplitude_step = amplitude / 2
last_amplitude = initial_amplitude * 2
last_ap = True

while abs(last_amplitude - amplitude) > epsilon or last_ap:
  print(f"Testing {amplitude}, current difference {round(abs(last_amplitude - amplitude), 3)} mV (epsilon = {epsilon} mV)")
  last_amplitude = amplitude

  simulation.ti.amplitude = amplitude
  simulation.run(duration = duration)

  last_ap = len(neuron.spike_times) > 0

  if last_ap:
    amplitude -= amplitude_step
  else:
    amplitude += amplitude_step

  amplitude_step = amplitude_step / 2

print(f"Sub-threshold amplitude for carrier {simulation.ti.carrier} with target {simulation.ti.beat}: {round(last_amplitude, 2)} V/m")

pd.DataFrame([{
  "carrier": simulation.ti.carrier,
  "beat": simulation.ti.beat,
  "phase": simulation.ti.phase,
  "phi": simulation.ti.phi_deg,
  "psi": simulation.ti.psi_deg,
  "amplitude": last_amplitude,
  "epsilon": epsilon,
  "amplitude_step": amplitude_step * 2,
  "last_ap": last_ap
}]).to_csv(snakemake.output["subthreshold_ef"])
