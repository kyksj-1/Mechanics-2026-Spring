from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ising_core import relaxation_trajectory, simulate_ising_observables

ROOT = Path(__file__).resolve().parent
ASSET = ROOT / "asset"
DOCS = ROOT / "docs"
ASSET.mkdir(exist_ok=True)
DOCS.mkdir(exist_ok=True)

TC = 2.0 / np.log(1.0 + np.sqrt(2.0))


def steady_energy_per_site(L: int, seed: int) -> float:
    e, _, _, _, n = simulate_ising_observables(
        L=L,
        T=TC,
        n_therm=12000,
        n_measure=25000,
        stride=2,
        seed=seed,
    )
    return e / n


def _single_run(args: tuple[int, int, int]) -> np.ndarray:
    L, steps, seed = args
    return relaxation_trajectory(L=L, T=TC, n_steps=steps, seed=seed)


def ensemble_relaxation(L: int, steps: int, n_runs: int, seed0: int) -> tuple[np.ndarray, np.ndarray]:
    jobs = [(L, steps, seed0 + i) for i in range(n_runs)]
    trajs = []
    with ProcessPoolExecutor() as ex:
        for tr in ex.map(_single_run, jobs):
            trajs.append(tr)
    arr = np.array(trajs)
    mean = np.mean(arr, axis=0)
    stderr = np.std(arr, axis=0, ddof=1) / np.sqrt(n_runs)
    return mean, stderr


def fit_power_law(ts: np.ndarray, ys: np.ndarray) -> tuple[float, float]:
    x = np.log(ts)
    y = np.log(ys)
    a, b = np.polyfit(x, y, 1)
    return a, b


def main() -> None:
    # L=16 energy relaxation
    steps16 = 1200
    mean_traj16, stderr_traj16 = ensemble_relaxation(L=16, steps=steps16, n_runs=220, seed0=41000)
    e_inf16 = steady_energy_per_site(16, seed=51000)
    e16_per_site = mean_traj16 / (16 * 16)
    e16_err_per_site = stderr_traj16 / (16 * 16)

    plt.figure(figsize=(7, 5))
    plt.plot(np.arange(steps16 + 1), e16_per_site, lw=1.4, label=r"$\langle E(t)\rangle/N$")
    plt.fill_between(
        np.arange(steps16 + 1),
        e16_per_site - e16_err_per_site,
        e16_per_site + e16_err_per_site,
        alpha=0.22,
        linewidth=0,
        label="standard error",
    )
    plt.axhline(e_inf16, color="r", ls="--", label=r"$\langle E(\infty)\rangle/N$")
    plt.xlabel("t (MCS)")
    plt.ylabel("Energy per site")
    plt.title("Relaxation at $L=16$, $T=T_c$")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ASSET / "b_relax_l16.png", dpi=180)
    plt.close()

    # finite-size effect
    sizes = [8, 16, 24, 32]
    steps = 1500
    runs = {8: 260, 16: 220, 24: 180, 32: 140}

    all_data = {}
    plt.figure(figsize=(7, 5))
    for L in sizes:
        mean_traj, stderr_traj = ensemble_relaxation(L=L, steps=steps, n_runs=runs[L], seed0=62000 + 1000 * L)
        e_inf = steady_energy_per_site(L, seed=73000 + L)
        delta = np.abs(mean_traj / (L * L) - e_inf)
        delta_err = stderr_traj / (L * L)
        t = np.arange(1, steps + 1)
        d = delta[1:]

        mask = (t >= 150) & (t <= 1200) & (d > 0)
        slope, intercept = fit_power_law(t[mask], d[mask])

        all_data[str(L)] = {
            "e_inf_per_site": float(e_inf),
            "slope": float(slope),
            "intercept": float(intercept),
            "delta_t": delta.tolist(),
            "delta_err_t": delta_err.tolist(),
        }

        plt.loglog(t, d, label=f"L={L}, slope={slope:.3f}")
        plt.fill_between(t, np.maximum(d - delta_err[1:], 1e-12), d + delta_err[1:], alpha=0.12, linewidth=0)

    plt.xlabel("t (MCS)")
    plt.ylabel(r"$\Delta(t)=|\langle E(t)\rangle/N-\langle E(\infty)\rangle/N|$")
    plt.title("Long-time relaxation at $T=T_c$")
    plt.grid(alpha=0.3, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(ASSET / "b_delta_loglog.png", dpi=180)
    plt.close()

    with (DOCS / "b_results.json").open("w", encoding="utf-8") as f:
        json.dump(
            {
                "Tc": float(TC),
                "L16": {
                    "e_inf_per_site": float(e_inf16),
                    "traj_per_site": e16_per_site.tolist(),
                    "traj_err_per_site": e16_err_per_site.tolist(),
                },
                "size_scan": all_data,
            },
            f,
            indent=2,
            ensure_ascii=False,
        )

    print(json.dumps({"Tc": float(TC), "L16_e_inf_per_site": float(e_inf16)}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
