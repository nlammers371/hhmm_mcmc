from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .gillespie import simulate_trace, TraceResult
from .model import PromoterModel, SimulationConfig


@dataclass(slots=True)
class DatasetResult:
    observed_fluo: np.ndarray
    true_states_sampled: np.ndarray
    traces: list[TraceResult]


def _sample_states_on_grid(transition_times: np.ndarray, states: np.ndarray, t_grid: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(transition_times, t_grid, side="right") - 1
    idx = np.clip(idx, 0, len(states) - 1)
    return states[idx]


def simulate_dataset(
    model: PromoterModel,
    config: SimulationConfig,
    n_traces: int,
    seed: int | None = None,
) -> DatasetResult:
    rng = np.random.default_rng(seed)
    observed = np.zeros((config.seq_length, n_traces), dtype=float)
    true_states = np.zeros((config.seq_length, n_traces), dtype=int)
    traces: list[TraceResult] = []

    t_grid = np.arange(config.seq_length, dtype=float) * config.delta_t
    for n in range(n_traces):
        tr = simulate_trace(model, config, rng=rng)
        traces.append(tr)
        observed[:, n] = tr.fluo_ms2
        true_states[:, n] = _sample_states_on_grid(tr.transition_times, tr.naive_states, t_grid)

    return DatasetResult(observed_fluo=observed, true_states_sampled=true_states, traces=traces)
