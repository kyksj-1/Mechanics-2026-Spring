from __future__ import annotations

import itertools
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ising_core import simulate_ising_measurement_series, simulate_ising_observables

ROOT = Path(__file__).resolve().parent
ASSET = ROOT / "asset"
DOCS = ROOT / "docs"
ASSET.mkdir(exist_ok=True)
DOCS.mkdir(exist_ok=True)
TC = 2.0 / np.log(1.0 + np.sqrt(2.0))


def exact_l4_t1() -> dict:
    L = 4
    N = L * L
    beta = 1.0
    Z = 0.0
    E_sum = 0.0
    states = 0

    for bits in itertools.product((-1, 1), repeat=N):
        s = np.array(bits, dtype=np.int8).reshape(L, L)
        e = 0
        for i in range(L):
            ip = (i + 1) % L
            for j in range(L):
                jp = (j + 1) % L
                e -= s[i, j] * s[ip, j]
                e -= s[i, j] * s[i, jp]
        w = np.exp(-beta * e)
        Z += w
        E_sum += e * w
        states += 1

    E = E_sum / Z
    F = -np.log(Z) / beta
    return {"L": L, "T": 1.0, "num_states": states, "E_exact": E, "F_exact": F, "E_per_site": E / N, "F_per_site": F / N}


def block_observables(
    energies: np.ndarray,
    magnetizations: np.ndarray,
    L: int,
    T: float,
    n_blocks: int = 24,
) -> dict:
    n_sites = L * L
    beta = 1.0 / T
    usable = (len(energies) // n_blocks) * n_blocks
    e_blocks = energies[:usable].reshape(n_blocks, -1)
    m_blocks = magnetizations[:usable].reshape(n_blocks, -1)

    rows = []
    for e, m in zip(e_blocks, m_blocks):
        e_mean = np.mean(e)
        e2_mean = np.mean(e * e)
        m2_mean = np.mean(m * m)
        absm_mean = np.mean(np.abs(m))
        rows.append(
            (
                e_mean / n_sites,
                m2_mean / (n_sites * n_sites),
                beta * beta * (e2_mean - e_mean * e_mean) / n_sites,
                beta * (m2_mean - absm_mean * absm_mean) / n_sites,
            )
        )

    arr = np.array(rows)
    means = np.mean(arr, axis=0)
    errs = np.std(arr, axis=0, ddof=1) / np.sqrt(n_blocks)
    return {
        "E_per_site": float(means[0]),
        "E_per_site_err": float(errs[0]),
        "m2": float(means[1]),
        "m2_err": float(errs[1]),
        "c": float(means[2]),
        "c_err": float(errs[2]),
        "chi": float(means[3]),
        "chi_err": float(errs[3]),
        "n_blocks": int(n_blocks),
        "block_size": int(usable // n_blocks),
    }


def _run_one_temperature(args: tuple[int, float, int]) -> dict:
    L, T, seed = args
    energies, magnetizations = simulate_ising_measurement_series(
        L=L,
        T=T,
        n_therm=4000,
        n_measure=12000,
        stride=3,
        seed=seed,
    )
    stats = block_observables(energies, magnetizations, L, T)
    stats.update({"L": int(L), "T": float(T)})
    return stats


def monte_carlo_scan() -> list[dict]:
    temps = np.round(np.arange(1.5, 3.0 + 1e-12, 0.1), 1)
    sizes = [8, 16, 32]
    jobs = []
    seed_base = 20260512
    for L in sizes:
        for k, T in enumerate(temps):
            jobs.append((L, float(T), seed_base + L * 1000 + k))

    results = []
    with ProcessPoolExecutor() as ex:
        for out in ex.map(_run_one_temperature, jobs):
            results.append(out)

    results.sort(key=lambda x: (x["L"], x["T"]))
    return results


def l4_mc_check() -> dict:
    energies, magnetizations = simulate_ising_measurement_series(
        L=4,
        T=1.0,
        n_therm=5000,
        n_measure=30000,
        stride=4,
        seed=2026,
    )
    stats = block_observables(energies, magnetizations, L=4, T=1.0, n_blocks=30)
    n_sites = 16
    return {
        "L": 4,
        "T": 1.0,
        "E_mc": stats["E_per_site"] * n_sites,
        "E_mc_err": stats["E_per_site_err"] * n_sites,
        "E_mc_per_site": stats["E_per_site"],
        "E_mc_per_site_err": stats["E_per_site_err"],
    }


def make_plots(results: list[dict]) -> None:
    fig_specs = [
        ("m2", "m2_err", r"$\langle m^2\rangle$", "a_m2_vs_t.png"),
        ("c", "c_err", r"$c$", "a_c_vs_t.png"),
        ("chi", "chi_err", r"$\chi$", "a_chi_vs_t.png"),
    ]

    for key, err_key, ylabel, fname in fig_specs:
        plt.figure(figsize=(7, 5))
        for L in [8, 16, 32]:
            sel = [r for r in results if r["L"] == L]
            temps = np.array([r["T"] for r in sel])
            vals = np.array([r[key] for r in sel])
            errs = np.array([r[err_key] for r in sel])
            plt.errorbar(temps, vals, yerr=errs, marker="o", capsize=3, lw=1.3, label=f"L={L}")
        plt.xlabel("T")
        plt.ylabel(ylabel)
        plt.title(ylabel + " vs T")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(ASSET / fname, dpi=180)
        plt.close()

    plt.figure(figsize=(7, 5))
    for L in [8, 16, 32]:
        sel = [r for r in results if r["L"] == L]
        temps = np.array([r["T"] for r in sel])
        x = (temps - TC) * L
        y = np.array([r["m2"] for r in sel]) * (L ** 0.25)
        yerr = np.array([r["m2_err"] for r in sel]) * (L ** 0.25)
        plt.errorbar(x, y, yerr=yerr, marker="o", capsize=3, lw=1.2, label=f"L={L}")
    plt.axvline(0.0, color="0.35", lw=1.0, ls=":")
    plt.xlabel(r"$(T-T_c)L^{1/\nu}$,  $\nu=1$")
    plt.ylabel(r"$L^{2\beta/\nu}\langle m^2\rangle$,  $2\beta/\nu=1/4$")
    plt.title("Finite-size data collapse of magnetization")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ASSET / "a_data_collapse_m2.png", dpi=180)
    plt.close()


def main() -> None:
    exact = exact_l4_t1()
    mc_l4 = l4_mc_check()
    scan = monte_carlo_scan()
    make_plots(scan)

    out = {
        "exact": exact,
        "mc_l4": mc_l4,
        "scan": [
            r
            for r in scan
        ],
    }

    with (DOCS / "a_results.json").open("w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print(json.dumps({"exact": exact, "mc_l4": mc_l4}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
