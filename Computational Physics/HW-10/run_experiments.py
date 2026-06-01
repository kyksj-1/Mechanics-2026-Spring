"""Run all numerical experiments for HW-10 and generate figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hopfield import (
    add_noise,
    cosine_similarity,
    energy_from_overlaps,
    overlaps,
    random_spins,
    zero_temperature_metropolis,
)


ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "asset"
DEFAULT_SEED = 20260601


def _round_list(values: np.ndarray, digits: int = 6) -> list[float]:
    return [round(float(value), digits) for value in values]


def _save_p1_figure(energy: np.ndarray, similarity: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)

    axes[0].plot(np.arange(len(energy)), energy, marker="o", color="#2a6f97")
    axes[0].set_xlabel("Sweep")
    axes[0].set_ylabel("Energy")
    axes[0].set_title("p=1 energy relaxation")
    axes[0].grid(alpha=0.25)

    axes[1].plot(np.arange(len(similarity)), similarity, marker="o", color="#c44536")
    axes[1].set_xlabel("Sweep")
    axes[1].set_ylabel("Similarity to stored pattern")
    axes[1].set_ylim(0.0, 1.05)
    axes[1].set_title("p=1 retrieval")
    axes[1].grid(alpha=0.25)

    fig.savefig(ASSET_DIR / "fig_p1_relaxation.png", dpi=180)
    plt.close(fig)


def _save_p10_figures(
    initial_similarity: np.ndarray,
    final_similarity: np.ndarray,
    final_overlap_matrix: np.ndarray,
    mean_energy: np.ndarray,
) -> None:
    index = np.arange(1, len(initial_similarity) + 1)

    fig, ax = plt.subplots(figsize=(10.5, 4.2), constrained_layout=True)
    width = 0.36
    ax.bar(index - width / 2, initial_similarity, width, label="Initial", color="#6a994e")
    ax.bar(index + width / 2, final_similarity, width, label="Final", color="#bc4749")
    ax.set_xlabel("Pattern index")
    ax.set_ylabel("Similarity to target pattern")
    ax.set_ylim(0.0, 1.08)
    ax.set_xticks(index)
    ax.set_title("N=2000, p=10 retrieval from delta=0.5 noise")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(ASSET_DIR / "fig_p10_similarity_bars.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 5.4), constrained_layout=True)
    image = ax.imshow(final_overlap_matrix, vmin=0.0, vmax=1.0, cmap="magma")
    ax.set_xlabel("Stored pattern index")
    ax.set_ylabel("Relaxed state index")
    ax.set_title("Absolute final overlap matrix")
    fig.colorbar(image, ax=ax, label="|overlap| / N")
    fig.savefig(ASSET_DIR / "fig_p10_overlap_matrix.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    ax.plot(np.arange(len(mean_energy)), mean_energy, marker="o", color="#386641")
    ax.set_xlabel("Sweep")
    ax.set_ylabel("Mean energy")
    ax.set_title("Mean energy over the 10 replicas")
    ax.grid(alpha=0.25)
    fig.savefig(ASSET_DIR / "fig_p10_mean_energy.png", dpi=180)
    plt.close(fig)


def _save_capacity_figures(summary: dict[str, dict[str, float]], samples: dict[int, np.ndarray]) -> None:
    p_values = sorted(samples)
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 6.0), constrained_layout=True)
    colors = ["#005f73", "#0a9396", "#94d2bd", "#ee9b00", "#ca6702", "#bb3e03", "#9b2226", "#5f0f40"]
    bins = np.linspace(0.0, 1.0, 21)

    for ax, color, p_value in zip(axes.flat, colors, p_values):
        values = samples[p_value]
        ax.hist(values, bins=bins, color=color, alpha=0.86, edgecolor="white")
        ax.axvline(np.mean(values), color="black", linewidth=1.3)
        ax.set_title(f"p={p_value}, mean={np.mean(values):.3f}")
        ax.set_xlim(0.0, 1.02)
        ax.set_xlabel("Final similarity")
        ax.set_ylabel("Count")
        ax.grid(axis="y", alpha=0.2)

    fig.savefig(ASSET_DIR / "fig_capacity_histograms.png", dpi=180)
    plt.close(fig)

    means = np.array([summary[str(p)]["mean"] for p in p_values])
    stds = np.array([summary[str(p)]["std"] for p in p_values])
    success = np.array([summary[str(p)]["success_rate_ge_0p9"] for p in p_values])

    fig, ax1 = plt.subplots(figsize=(8.2, 4.4), constrained_layout=True)
    ax1.errorbar(
        p_values,
        means,
        yerr=stds,
        marker="o",
        linewidth=2.0,
        capsize=4,
        color="#1d3557",
        label="Mean similarity",
    )
    ax1.set_xlabel("Number of stored patterns p")
    ax1.set_ylabel("Final similarity")
    ax1.set_ylim(0.0, 1.05)
    ax1.grid(alpha=0.25)

    ax2 = ax1.twinx()
    ax2.plot(p_values, success, marker="s", linewidth=2.0, color="#e76f51", label="S >= 0.9")
    ax2.set_ylabel("Retrieval success rate")
    ax2.set_ylim(0.0, 1.05)

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, frameon=False, loc="lower left")
    ax1.set_title("Capacity trend for N=2000")
    fig.savefig(ASSET_DIR / "fig_capacity_summary.png", dpi=180)
    plt.close(fig)


def experiment_p1(rng: np.random.Generator) -> dict[str, object]:
    n_sites = 1000
    patterns = random_spins(rng, (1, n_sites))
    initial_state = random_spins(rng, n_sites)

    initial_overlap = overlaps(patterns, initial_state)
    initial_similarity = cosine_similarity(initial_state, patterns[0])[0]
    initial_energy = energy_from_overlaps(initial_overlap, n_sites, 1)[0]

    result = zero_temperature_metropolis(
        patterns,
        initial_state,
        rng,
        max_sweeps=20,
        record_against=patterns,
    )
    final_state = result.states[0]
    final_similarity = cosine_similarity(final_state, patterns[0])[0]
    final_energy = energy_from_overlaps(result.overlaps, n_sites, 1)[0]
    if np.array_equal(final_state, patterns[0]):
        final_state_type = "+v"
    elif np.array_equal(final_state, -patterns[0]):
        final_state_type = "-v"
    else:
        final_state_type = "not exactly +/-v"

    _save_p1_figure(result.energy_history[:, 0], result.similarity_history[:, 0])

    return {
        "N": n_sites,
        "p": 1,
        "initial_similarity": float(initial_similarity),
        "initial_energy": float(initial_energy),
        "final_similarity": float(final_similarity),
        "final_energy": float(final_energy),
        "ground_energy_exact": float(-(n_sites - 1) / 2),
        "sweeps": result.sweeps,
        "flips_per_sweep": result.flips_per_sweep,
        "final_state_type": final_state_type,
    }


def experiment_p10(rng: np.random.Generator) -> dict[str, object]:
    n_sites = 2000
    n_patterns = 10
    delta = 0.5
    patterns = random_spins(rng, (n_patterns, n_sites))
    initial_states = add_noise(rng, patterns, delta)

    initial_similarity = cosine_similarity(initial_states, patterns)
    result = zero_temperature_metropolis(
        patterns,
        initial_states,
        rng,
        max_sweeps=25,
        record_against=patterns,
    )
    final_similarity = cosine_similarity(result.states, patterns)
    final_overlap_matrix = np.abs(result.overlaps) / n_sites
    final_energy = energy_from_overlaps(result.overlaps, n_sites, n_patterns)

    _save_p10_figures(
        initial_similarity,
        final_similarity,
        final_overlap_matrix,
        result.energy_history.mean(axis=1),
    )

    return {
        "N": n_sites,
        "p": n_patterns,
        "delta": delta,
        "initial_similarity": _round_list(initial_similarity),
        "final_similarity": _round_list(final_similarity),
        "final_energy": _round_list(final_energy),
        "mean_final_similarity": float(np.mean(final_similarity)),
        "min_final_similarity": float(np.min(final_similarity)),
        "max_off_diagonal_abs_overlap": float(
            np.max(final_overlap_matrix - np.diag(np.diag(final_overlap_matrix)))
        ),
        "sweeps": result.sweeps,
        "flips_per_sweep": result.flips_per_sweep,
    }


def experiment_capacity(
    rng: np.random.Generator,
    p_values: list[int],
    repeats: int,
    max_sweeps: int,
) -> tuple[dict[str, dict[str, float]], dict[int, np.ndarray]]:
    n_sites = 2000
    delta = 0.5
    samples: dict[int, np.ndarray] = {}
    summary: dict[str, dict[str, float]] = {}

    for p_value in p_values:
        p_samples: list[np.ndarray] = []
        sweep_counts: list[int] = []
        last_flips: list[int] = []
        for _ in range(repeats):
            patterns = random_spins(rng, (p_value, n_sites))
            initial_states = add_noise(rng, patterns, delta)
            result = zero_temperature_metropolis(
                patterns,
                initial_states,
                rng,
                max_sweeps=max_sweeps,
                record_against=None,
            )
            final_similarity = cosine_similarity(result.states, patterns)
            p_samples.append(final_similarity)
            sweep_counts.append(result.sweeps)
            last_flips.append(result.flips_per_sweep[-1])

        values = np.concatenate(p_samples)
        samples[p_value] = values
        summary[str(p_value)] = {
            "samples": int(values.size),
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "median": float(np.median(values)),
            "q10": float(np.quantile(values, 0.10)),
            "q90": float(np.quantile(values, 0.90)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "success_rate_ge_0p9": float(np.mean(values >= 0.9)),
            "mean_sweeps": float(np.mean(sweep_counts)),
            "max_last_sweep_flips": int(np.max(last_flips)),
        }

    _save_capacity_figures(summary, samples)
    return summary, samples


def print_results(results: dict[str, object]) -> None:
    p1 = results["p1"]
    p10 = results["p10"]
    capacity = results["capacity"]

    print(f"SEED {results['seed']}")
    print(f"RUNTIME_SECONDS {results['runtime_seconds']:.3f}")
    print("P1")
    print(f"N {p1['N']} p {p1['p']}")
    print(f"initial_similarity {p1['initial_similarity']:.6f}")
    print(f"initial_energy {p1['initial_energy']:.6f}")
    print(f"final_similarity {p1['final_similarity']:.6f}")
    print(f"final_energy {p1['final_energy']:.6f}")
    print(f"ground_energy_exact {p1['ground_energy_exact']:.6f}")
    print(f"sweeps {p1['sweeps']}")
    print(f"flips_per_sweep {p1['flips_per_sweep']}")
    print(f"final_state_type {p1['final_state_type']}")

    print("P10")
    print(f"N {p10['N']} p {p10['p']} delta {p10['delta']}")
    print(f"initial_similarity {p10['initial_similarity']}")
    print(f"final_similarity {p10['final_similarity']}")
    print(f"final_energy {p10['final_energy']}")
    print(f"mean_final_similarity {p10['mean_final_similarity']:.6f}")
    print(f"min_final_similarity {p10['min_final_similarity']:.6f}")
    print(f"max_off_diagonal_abs_overlap {p10['max_off_diagonal_abs_overlap']:.6f}")
    print(f"sweeps {p10['sweeps']}")
    print(f"flips_per_sweep {p10['flips_per_sweep']}")

    print("CAPACITY")
    print("p samples mean std median q10 q90 min max success_ge_0p9 mean_sweeps max_last_sweep_flips")
    for p_key in sorted(capacity, key=lambda item: int(item)):
        row = capacity[p_key]
        print(
            f"{p_key} {row['samples']} {row['mean']:.6f} {row['std']:.6f} "
            f"{row['median']:.6f} {row['q10']:.6f} {row['q90']:.6f} "
            f"{row['min']:.6f} {row['max']:.6f} "
            f"{row['success_rate_ge_0p9']:.6f} {row['mean_sweeps']:.3f} "
            f"{row['max_last_sweep_flips']}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--max-sweeps", type=int, default=100)
    parser.add_argument(
        "--p-values",
        type=int,
        nargs="+",
        default=[50, 100, 150, 200, 250, 300, 350, 400],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ASSET_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(args.seed)
    start = perf_counter()

    p1 = experiment_p1(rng)
    p10 = experiment_p10(rng)
    capacity, _ = experiment_capacity(
        rng,
        p_values=args.p_values,
        repeats=args.repeats,
        max_sweeps=args.max_sweeps,
    )

    results: dict[str, object] = {
        "seed": args.seed,
        "runtime_seconds": perf_counter() - start,
        "p1": p1,
        "p10": p10,
        "capacity": capacity,
    }

    with (ASSET_DIR / "numeric_results.json").open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print_results(results)


if __name__ == "__main__":
    main()
