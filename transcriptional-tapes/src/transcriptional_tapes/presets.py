from __future__ import annotations

import numpy as np

from .model import PromoterModel


def stationary_pi_from_rate_matrix(rate_matrix: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(rate_matrix)
    idx = int(np.argmin(np.abs(vals)))
    pi = np.real(vecs[:, idx])
    pi = pi / np.sum(pi)
    if np.any(pi < 0):
        pi = np.abs(pi)
        pi = pi / np.sum(pi)
    return pi


def two_state_default() -> PromoterModel:
    r = np.array([[-0.015, 0.02], [0.015, -0.02]], dtype=float)
    v = np.array([0.0, 4.0], dtype=float) / 20.0
    pi0 = stationary_pi_from_rate_matrix(r)
    return PromoterModel(rate_matrix=r, emission_rates=v, sigma=1.5, pi0=pi0)


def three_state_default() -> PromoterModel:
    r = np.array(
        [[-0.0200, 0.0210, 0.0], [0.0200, -0.0330, 0.0700], [0.0, 0.0120, -0.0700]],
        dtype=float,
    )
    v = np.array([0.0, 2.0, 5.0], dtype=float) / 20.0
    pi0 = stationary_pi_from_rate_matrix(r)
    return PromoterModel(rate_matrix=r, emission_rates=v, sigma=1.5, pi0=pi0)
