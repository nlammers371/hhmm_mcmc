from __future__ import annotations


def ms2_loading_coeff_integral(alpha: float, w: int, delta_t: float, t1: float, t2: float) -> float:
    """Integral of MS2 loading coefficient over [t1, t2].

    Port of the piecewise form used in MATLAB `ms2_loading_coeff_integral`.
    """
    eps = 1e-8
    if t1 < -eps or t2 < -eps or (t1 - t2) > eps or (t2 - w * delta_t) > eps:
        raise ValueError("require 0 <= t1 <= t2 <= w*delta_t")
    if alpha < 0 or alpha - w > eps:
        raise ValueError("require 0 <= alpha <= w")

    t_ms2 = alpha * delta_t

    if t2 <= t_ms2:
        return (t2 * t2 - t1 * t1) / (2.0 * t_ms2) if t_ms2 > 0 else 0.0
    if t1 >= t_ms2:
        return t2 - t1

    coeff = (t_ms2 * t_ms2 - t1 * t1) / (2.0 * t_ms2) if t_ms2 > 0 else 0.0
    coeff += t2 - t_ms2
    return coeff
