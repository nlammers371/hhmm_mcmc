import numpy as np

from transcriptional_tapes.batch import simulate_dataset
from transcriptional_tapes.model import SimulationConfig
from transcriptional_tapes.ms2 import ms2_loading_coeff_integral
from transcriptional_tapes.presets import two_state_default


def test_ms2_integral_regions():
    alpha, w, dt = 2.0, 7, 1.0
    # entirely in ramp
    v1 = ms2_loading_coeff_integral(alpha, w, dt, 0.2, 1.2)
    assert np.isclose(v1, (1.2**2 - 0.2**2) / (2 * alpha * dt))
    # entirely post-ramp
    v2 = ms2_loading_coeff_integral(alpha, w, dt, 2.2, 4.2)
    assert np.isclose(v2, 2.0)


def test_dataset_shapes_and_reproducibility():
    model = two_state_default()
    config = SimulationConfig(seq_length=50, delta_t=20.0, memory_steps=7, alpha=1.5)
    out1 = simulate_dataset(model, config, n_traces=5, seed=123)
    out2 = simulate_dataset(model, config, n_traces=5, seed=123)

    assert out1.observed_fluo.shape == (50, 5)
    assert out1.true_states_sampled.shape == (50, 5)
    assert np.allclose(out1.observed_fluo, out2.observed_fluo)
