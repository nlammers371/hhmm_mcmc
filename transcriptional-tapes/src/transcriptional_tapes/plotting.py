from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_single_trace(time_axis: np.ndarray, fluo_ms2: np.ndarray, fluo_raw: np.ndarray | None = None):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(time_axis, fluo_ms2, label="MS2 fluorescence", lw=2)
    if fluo_raw is not None:
        ax.plot(time_axis, fluo_raw, label="raw fluorescence", alpha=0.7)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("fluorescence (a.u.)")
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig, ax


def plot_ensemble_mean(time_axis: np.ndarray, trajectories: np.ndarray):
    mean = trajectories.mean(axis=1)
    std = trajectories.std(axis=1)
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(time_axis, mean, lw=2, color="tab:blue")
    ax.fill_between(time_axis, mean - std, mean + std, color="tab:blue", alpha=0.2)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("fluorescence (a.u.)")
    ax.set_title("Ensemble mean ± 1 SD")
    fig.tight_layout()
    return fig, ax
