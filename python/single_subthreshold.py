

from python.ca1 import NeuronCA1, Simulation

neuron = NeuronCA1()
simulation = Simulation(neurons = [ neuron ])

simulation.ti.beat = 5 # Hz - TODO: take from snakemake
simulation.ti.carrier = 1000 # Hz - TODO: take from snakemake
simulation.ti.phase = 10
simulation.ti.delay = 0

simulation.ti.phi_deg = 90
simulation.ti.psi_deg = 90

# Get the optimal stimulation duration to include a minimum of 10 5Hz oscillations
duration = simulation.ti.get_optimal_duration(min_oscillations = 10)
simulation.ti.duration = duration

initial_amplitude = 2000 # V/m
epsilon = 0.5 # V/m

amplitude = initial_amplitude
amplitude_step = amplitude / 2
last_amplitude = initial_amplitude * 2
last_ap = True

while abs(last_amplitude - amplitude) > epsilon or last_ap:
    print(f"Testing {amplitude}")
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
