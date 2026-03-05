from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(slots=True)
class PromoterModel:
    """Continuous-time promoter model and fluorescence parameters."""

    rate_matrix: np.ndarray
    emission_rates: np.ndarray
    sigma: float
    pi0: np.ndarray

    def validate(self) -> None:
        k = self.rate_matrix.shape[0]
        if self.rate_matrix.shape != (k, k):
            raise ValueError("rate_matrix must be square")
        if self.emission_rates.shape != (k,):
            raise ValueError("emission_rates must have shape (K,)")
        if self.pi0.shape != (k,):
            raise ValueError("pi0 must have shape (K,)")
        if not np.isclose(np.sum(self.pi0), 1.0):
            raise ValueError("pi0 must sum to 1")
        if np.any(self.pi0 < 0):
            raise ValueError("pi0 must be non-negative")
        if np.any(np.diag(self.rate_matrix) >= 0):
            raise ValueError("diagonal rates must be negative")
        if np.any(self.rate_matrix - np.diag(np.diag(self.rate_matrix)) < 0):
            raise ValueError("off-diagonal rates must be non-negative")


@dataclass(slots=True)
class SimulationConfig:
    """Configuration for uniformly sampled fluorescence trajectories."""

    seq_length: int
    delta_t: float
    memory_steps: int
    alpha: float

    def validate(self) -> None:
        if self.seq_length <= 0:
            raise ValueError("seq_length must be > 0")
        if self.delta_t <= 0:
            raise ValueError("delta_t must be > 0")
        if self.memory_steps <= 0:
            raise ValueError("memory_steps must be > 0")
        if not 0 <= self.alpha <= self.memory_steps:
            raise ValueError("alpha must satisfy 0 <= alpha <= memory_steps")
