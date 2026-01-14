# Simulating the effects of temporal interference stimulation on a two-neuron network

## How to run the simulation

1. Create and activate `neuron` conda environment.

2. Compile NMODL files:

```{shell}
nrnivmodl nrnmod
```

3. Run snakemake

## Tasks

### Jan 2026

**Step 1.** Create a two-neuron network.
- [ ] Find out how to duplicate and translate the morphology.
- [ ] Make sure both neurons react to external electric fields.
- [ ] Set up a method to comfortably translate the morphology (shift to side, rotate on both axes).
- [ ] Connect both neurons through the soma (unidirectional).

**Step 2.** Set up the parameters of the two-neuron network.
- [ ] If neuron 1 connects to neuron 2, neuron 2 should spike right after neuron 1 spikes (threshold should be 0, weight should be high enough to generate a spike).
- [ ] Test unconnected network: if both neurons behave the same, rotate one of them. They should react to TI differently due to different position in the electric field.

**Step 3.** Evaluate subthreshold electric field amplitude for different parameter settings.
- [ ] Test 1kHz and 9kHz carriers for 5Hz beat frequency. Evalute subthreshold amplitude dependence on the distance between the neurons and their relative position.
- [ ] Subthreshold amplitude is defined as the minimum electric field amplitude when both neurons spike in a time window of 3 ms.
- [ ] Test unidirectional connection in both directions.
