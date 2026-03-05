# transcriptional-tapes

Python-first simulation toolkit for generating synthetic transcriptional fluorescence trajectories from continuous-time promoter dynamics with MS2 signal convolution.

## Conceptual overview (Gillespie simulation assumptions)

This package ports the simulation side of the MATLAB workflow and keeps the same modeling assumptions:

1. **Promoter switching is a continuous-time Markov chain (CTMC)** over `K` promoter states with rate matrix `R` (sec⁻¹).
2. **Event times are sampled by Gillespie SSA** using exponential waiting times from the current state's exit rate.
3. **State-dependent emission rates** convert promoter states into polymerase initiation intensity.
4. **MS2 fluorescence has memory**: observed signal depends on recent promoter activity over a finite window `w`.
5. **Finite MS2 loop loading is explicit** via a piecewise integral kernel controlled by `alpha` (fractional loading period).
6. **Sampling is discrete in acquisition time** (`delta_t`) even though latent state dynamics are continuous-time.
7. **Observation noise is additive Gaussian**, producing noisy traces from noiseless trajectories.
8. **Batch simulation returns both observables and latent summaries**, enabling downstream benchmarking and plotting.

## Install

```bash
pip install -e .
```

## Quick usage

```python
import numpy as np
from transcriptional_tapes.model import SimulationConfig
from transcriptional_tapes.presets import three_state_default
from transcriptional_tapes.batch import simulate_dataset

model = three_state_default()
config = SimulationConfig(seq_length=120, delta_t=20.0, memory_steps=7, alpha=30/140 * 7)

out = simulate_dataset(model=model, config=config, n_traces=50, seed=7)
print(out.observed_fluo.shape)       # (120, 50)
print(out.true_states_sampled.shape) # (120, 50)
```

## Included modules

- `transcriptional_tapes.gillespie`: single-trace simulation from CTMC + MS2 kernel.
- `transcriptional_tapes.batch`: helper for generating trajectory ensembles.
- `transcriptional_tapes.plotting`: lightweight plotting helpers.
- `transcriptional_tapes.presets`: 2-state and 3-state default models inspired by the original MATLAB settings.

## Notebook

See `notebooks/simulation_quickstart.ipynb` for an end-to-end example.
