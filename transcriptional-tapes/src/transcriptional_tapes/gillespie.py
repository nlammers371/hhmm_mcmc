from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .model import PromoterModel, SimulationConfig
from .ms2 import ms2_loading_coeff_integral


@dataclass(slots=True)
class TraceResult:
    fluo: np.ndarray
    fluo_ms2: np.ndarray
    fluo_no_noise: np.ndarray
    fluo_ms2_no_noise: np.ndarray
    transition_times: np.ndarray
    naive_states: np.ndarray


def simulate_trace(
    model: PromoterModel,
    config: SimulationConfig,
    rng: np.random.Generator | None = None,
) -> TraceResult:
    """Simulate one fluorescence trajectory from CTMC promoter dynamics."""
    model.validate()
    config.validate()
    if rng is None:
        rng = np.random.default_rng()

    k = model.rate_matrix.shape[0]
    t_max = config.seq_length * config.delta_t
    times_unif = np.arange(1, config.seq_length + 1, dtype=float) * config.delta_t

    state = rng.choice(np.arange(k), p=model.pi0)
    states = [state]
    times = [0.0]
    t = 0.0

    while t < t_max:
        lam = -model.rate_matrix[state, state]
        dt = rng.exponential(1.0 / lam)
        t += dt
        rates = model.rate_matrix[:, state].copy()
        rates[state] = 0.0
        probs = rates / lam
        state = int(rng.choice(np.arange(k), p=probs))
        states.append(state)
        times.append(t)

    transition_times = np.asarray(times, dtype=float)
    naive_states = np.asarray(states, dtype=int)

    fluo = np.zeros(config.seq_length, dtype=float)
    fluo_ms2 = np.zeros(config.seq_length, dtype=float)

    for idx, t_end in enumerate(times_unif):
        t_start = max(0.0, t_end - config.memory_steps * config.delta_t)
        i_start = np.searchsorted(transition_times, t_start, side="right") - 1
        i_end = np.searchsorted(transition_times, t_end, side="right") - 1

        times_window = transition_times[i_start : i_end + 2].copy()
        times_window[0] = t_start
        times_window[-1] = t_end
        naive_window = naive_states[i_start : i_end + 2]

        for i in range(len(naive_window) - 1):
            t1 = t_end - times_window[i + 1]
            t2 = t_end - times_window[i]
            emit = model.emission_rates[naive_window[i]]
            fluo_ms2[idx] += emit * ms2_loading_coeff_integral(
                config.alpha, config.memory_steps, config.delta_t, t1, t2
            )
            fluo[idx] += emit * (t2 - t1)

    noise = rng.normal(0.0, model.sigma, size=config.seq_length)
    return TraceResult(
        fluo=fluo + noise,
        fluo_ms2=fluo_ms2 + noise,
        fluo_no_noise=fluo,
        fluo_ms2_no_noise=fluo_ms2,
        transition_times=transition_times,
        naive_states=naive_states,
    )
