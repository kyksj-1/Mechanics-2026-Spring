from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numba import njit


@dataclass(frozen=True)
class SimParams:
    L: int
    T: float
    n_therm: int
    n_measure: int
    stride: int
    seed: int


@njit(cache=True)
def random_spins(L: int, rng_seed: int) -> np.ndarray:
    np.random.seed(rng_seed)
    spins = np.random.randint(0, 2, size=(L, L)).astype(np.int8)
    spins = 2 * spins - 1
    return spins


@njit(cache=True)
def total_energy(spins: np.ndarray) -> int:
    L = spins.shape[0]
    e = 0
    for i in range(L):
        ip = (i + 1) % L
        for j in range(L):
            jp = (j + 1) % L
            s = spins[i, j]
            e -= s * spins[ip, j]
            e -= s * spins[i, jp]
    return e


@njit(cache=True)
def total_magnetization(spins: np.ndarray) -> int:
    return int(np.sum(spins))


@njit(cache=True)
def _sweep_metropolis(
    spins: np.ndarray,
    L: int,
    beta: float,
    exp_lut: np.ndarray,
    energy: int,
    magnetization: int,
) -> Tuple[int, int]:
    n_sites = L * L
    for _ in range(n_sites):
        i = np.random.randint(0, L)
        j = np.random.randint(0, L)
        s = spins[i, j]
        nn = (
            spins[(i + 1) % L, j]
            + spins[(i - 1) % L, j]
            + spins[i, (j + 1) % L]
            + spins[i, (j - 1) % L]
        )
        dE = 2 * s * nn
        if dE <= 0:
            spins[i, j] = -s
            energy += dE
            magnetization += -2 * s
        else:
            idx = dE // 4
            if np.random.random() < exp_lut[idx]:
                spins[i, j] = -s
                energy += dE
                magnetization += -2 * s
    return energy, magnetization


@njit(cache=True)
def simulate_ising_observables(
    L: int,
    T: float,
    n_therm: int,
    n_measure: int,
    stride: int,
    seed: int,
) -> Tuple[float, float, float, float, float]:
    beta = 1.0 / T
    spins = random_spins(L, seed)
    energy = total_energy(spins)
    magnetization = total_magnetization(spins)

    exp_lut = np.zeros(3, dtype=np.float64)
    exp_lut[1] = math.exp(-beta * 4.0)
    exp_lut[2] = math.exp(-beta * 8.0)

    for _ in range(n_therm):
        energy, magnetization = _sweep_metropolis(spins, L, beta, exp_lut, energy, magnetization)

    e_sum = 0.0
    e2_sum = 0.0
    m2_sum = 0.0
    abs_m_sum = 0.0

    for _ in range(n_measure):
        for _ in range(stride):
            energy, magnetization = _sweep_metropolis(spins, L, beta, exp_lut, energy, magnetization)

        e = float(energy)
        m = float(magnetization)
        e_sum += e
        e2_sum += e * e
        m2_sum += m * m
        abs_m_sum += abs(m)

    norm = 1.0 / n_measure
    return e_sum * norm, e2_sum * norm, m2_sum * norm, abs_m_sum * norm, float(L * L)


@njit(cache=True)
def simulate_ising_measurement_series(
    L: int,
    T: float,
    n_therm: int,
    n_measure: int,
    stride: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    beta = 1.0 / T
    spins = random_spins(L, seed)
    energy = total_energy(spins)
    magnetization = total_magnetization(spins)

    exp_lut = np.zeros(3, dtype=np.float64)
    exp_lut[1] = math.exp(-beta * 4.0)
    exp_lut[2] = math.exp(-beta * 8.0)

    for _ in range(n_therm):
        energy, magnetization = _sweep_metropolis(spins, L, beta, exp_lut, energy, magnetization)

    energies = np.zeros(n_measure, dtype=np.float64)
    magnetizations = np.zeros(n_measure, dtype=np.float64)
    for k in range(n_measure):
        for _ in range(stride):
            energy, magnetization = _sweep_metropolis(spins, L, beta, exp_lut, energy, magnetization)
        energies[k] = energy
        magnetizations[k] = magnetization

    return energies, magnetizations


@njit(cache=True)
def relaxation_trajectory(
    L: int,
    T: float,
    n_steps: int,
    seed: int,
) -> np.ndarray:
    beta = 1.0 / T
    spins = random_spins(L, seed)
    energy = total_energy(spins)
    magnetization = total_magnetization(spins)

    exp_lut = np.zeros(3, dtype=np.float64)
    exp_lut[1] = math.exp(-beta * 4.0)
    exp_lut[2] = math.exp(-beta * 8.0)

    traj = np.zeros(n_steps + 1, dtype=np.float64)
    traj[0] = energy
    for t in range(1, n_steps + 1):
        energy, magnetization = _sweep_metropolis(spins, L, beta, exp_lut, energy, magnetization)
        traj[t] = energy
    return traj
