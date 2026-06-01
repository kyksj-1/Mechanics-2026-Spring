"""Utilities for the fully connected Ising/Hopfield model in HW-10.

The implementation keeps the pattern overlaps

    a_mu = sum_i v_i^mu s_i

as the dynamical variables.  This avoids constructing the dense coupling
matrix and makes each single-spin update O(p) instead of O(N).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


SpinArray = np.ndarray


@dataclass
class EvolutionResult:
    """Container returned by the zero-temperature dynamics."""

    states: SpinArray
    overlaps: np.ndarray
    sweeps: int
    flips_per_sweep: list[int]
    energy_history: np.ndarray
    similarity_history: np.ndarray


def random_spins(
    rng: np.random.Generator,
    shape: int | tuple[int, ...],
    dtype: np.dtype = np.int8,
) -> SpinArray:
    """Draw independent Ising spins with values -1 and +1."""

    return rng.choice(np.array([-1, 1], dtype=dtype), size=shape)


def add_noise(
    rng: np.random.Generator,
    patterns: SpinArray,
    delta: float,
) -> SpinArray:
    """Return noisy copies of the input patterns.

    At each site, with probability ``delta`` the spin is replaced by a new
    random Ising value.  This matches the pseudocode in the problem set.
    """

    noisy = patterns.copy()
    mask = rng.random(patterns.shape) < delta
    noisy[mask] = random_spins(rng, int(mask.sum()))
    return noisy.astype(np.int8, copy=False)


def overlaps(patterns: SpinArray, states: SpinArray) -> np.ndarray:
    """Compute overlap matrix A[r, mu] = sum_i s_i^r v_i^mu."""

    states_2d = np.atleast_2d(states).astype(np.float32, copy=False)
    patterns_2d = np.atleast_2d(patterns).astype(np.float32, copy=False)
    return states_2d @ patterns_2d.T


def energy_from_overlaps(overlap_values: np.ndarray, n_sites: int, n_patterns: int) -> np.ndarray:
    """Energy from overlaps for J_ii = 0.

    E = -1/(2N) * sum_mu [(sum_i v_i^mu s_i)^2 - N]
      = -1/(2N) * sum_mu a_mu^2 + p/2.
    """

    overlap_2d = np.atleast_2d(overlap_values)
    return -0.5 * np.sum(overlap_2d * overlap_2d, axis=1) / n_sites + 0.5 * n_patterns


def cosine_similarity(states: SpinArray, patterns: SpinArray) -> np.ndarray:
    """Absolute cosine similarity between rows of ``states`` and ``patterns``."""

    state_rows = np.atleast_2d(states).astype(np.float32, copy=False)
    pattern_rows = np.atleast_2d(patterns).astype(np.float32, copy=False)
    if state_rows.shape != pattern_rows.shape:
        raise ValueError("states and patterns must have matching row-wise shapes")
    return np.abs(np.sum(state_rows * pattern_rows, axis=1)) / state_rows.shape[1]


def zero_temperature_metropolis(
    patterns: SpinArray,
    initial_states: SpinArray,
    rng: np.random.Generator,
    max_sweeps: int = 30,
    record_against: SpinArray | None = None,
) -> EvolutionResult:
    """Run asynchronous beta -> infinity single-spin Metropolis dynamics.

    The update is exact single-spin dynamics.  Independent replicas are evolved
    in the same random site order and vectorized over the replica index.
    A proposed flip is accepted only when Delta E < 0; Delta E = 0 moves do not
    change the energy and are omitted to make convergence well defined.
    """

    pattern_rows = np.atleast_2d(patterns).astype(np.int8, copy=False)
    states = np.atleast_2d(initial_states).astype(np.int8, copy=True)
    n_patterns, n_sites = pattern_rows.shape
    pattern_float = pattern_rows.astype(np.float32)
    overlap_values = overlaps(pattern_rows, states)

    energy_history: list[np.ndarray] = [
        energy_from_overlaps(overlap_values, n_sites, n_patterns)
    ]
    similarity_history: list[np.ndarray] = []
    if record_against is not None:
        similarity_history.append(cosine_similarity(states, record_against))

    flips_per_sweep: list[int] = []
    for sweep in range(max_sweeps):
        flips_this_sweep = 0
        for site in rng.permutation(n_sites):
            column = pattern_float[:, site]
            local_field_numerator = overlap_values @ column - n_patterns * states[:, site]
            delta_energy_numerator = 2.0 * states[:, site] * local_field_numerator
            flip_mask = delta_energy_numerator < 0.0
            if not np.any(flip_mask):
                continue

            old_spins = states[flip_mask, site].astype(np.float32)
            states[flip_mask, site] *= -1
            overlap_values[flip_mask] -= 2.0 * old_spins[:, None] * column[None, :]
            flips_this_sweep += int(np.count_nonzero(flip_mask))

        flips_per_sweep.append(flips_this_sweep)
        energy_history.append(energy_from_overlaps(overlap_values, n_sites, n_patterns))
        if record_against is not None:
            similarity_history.append(cosine_similarity(states, record_against))
        if flips_this_sweep == 0:
            return EvolutionResult(
                states=states,
                overlaps=overlap_values,
                sweeps=sweep + 1,
                flips_per_sweep=flips_per_sweep,
                energy_history=np.vstack(energy_history),
                similarity_history=np.vstack(similarity_history)
                if similarity_history
                else np.empty((0, states.shape[0])),
            )

    return EvolutionResult(
        states=states,
        overlaps=overlap_values,
        sweeps=max_sweeps,
        flips_per_sweep=flips_per_sweep,
        energy_history=np.vstack(energy_history),
        similarity_history=np.vstack(similarity_history)
        if similarity_history
        else np.empty((0, states.shape[0])),
    )
