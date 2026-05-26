"""
Part B.1: L=12 TFIM 基态/第一激发态能量与 <sigma_1^x>

对 h = 0.5, 1.0, 2.0 三个横场强度分别对角化, 输出能量和磁化期望.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh

from tfim_core import build_tfim_hamiltonian, sigma_x_op, expectation


ASSET_DIR = os.path.join(os.path.dirname(__file__), "asset")
os.makedirs(ASSET_DIR, exist_ok=True)


def run_b1():
    L = 12
    J = 1.0
    h_list = [0.5, 1.0, 2.0]

    print("=" * 70)
    print(f"Part B.1: L = {L}, J = {J} 横场 Ising 模型基态/第一激发态")
    print("=" * 70)
    print(f"{'h':>6} | {'E_0':>14} | {'E_1':>14} | {'gap':>12} | {'<sigma_1^x>_GS':>16}")
    print("-" * 72)

    results = {}
    Sx0 = sigma_x_op(L, site=0)  # 第一个格点 (i=1 in 1-indexed)

    for h in h_list:
        H = build_tfim_hamiltonian(L, J, h, pbc=True)
        # 求最低 2 个本征值/向量
        evals, evecs = eigsh(H, k=2, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]
        E0, E1 = evals[0], evals[1]
        psi0 = evecs[:, 0].astype(np.complex128)
        sx_exp = expectation(psi0, Sx0).real
        results[h] = (E0, E1, sx_exp)
        print(f"{h:>6.2f} | {E0:>14.8f} | {E1:>14.8f} | {E1 - E0:>12.8f} | {sx_exp:>16.8f}")

    # 可视化: 三个 h 下基态在 Z 基的概率分布 (前 64 个最大)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, h in zip(axes, h_list):
        H = build_tfim_hamiltonian(L, J, h, pbc=True)
        evals, evecs = eigsh(H, k=1, which='SA')
        prob = np.abs(evecs[:, 0]) ** 2
        idx_sorted = np.argsort(prob)[::-1][:32]
        ax.bar(range(len(idx_sorted)), prob[idx_sorted], color='steelblue', edgecolor='k', linewidth=0.4)
        ax.set_title(f"h = {h}    E_0 = {evals[0]:.4f}")
        ax.set_xlabel("basis index (sorted by prob)")
        ax.set_ylabel("|c|^2")
        ax.set_yscale('log')
        ax.grid(alpha=0.3)
    fig.suptitle(f"TFIM ground-state amplitudes in Z basis (L={L})")
    fig.tight_layout()
    out = os.path.join(ASSET_DIR, "b1_ground_state_amplitudes.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\n[Saved] {out}")

    # 可视化: <sigma_1^x>_GS 随 h 扫描, 看相变
    print("\n扫描 h ∈ [0.1, 3.0] 观察相变行为...")
    h_scan = np.linspace(0.1, 3.0, 30)
    sx_scan = np.empty_like(h_scan)
    E0_scan = np.empty_like(h_scan)
    gap_scan = np.empty_like(h_scan)
    for k, hh in enumerate(h_scan):
        H = build_tfim_hamiltonian(L, J, hh, pbc=True)
        evals, evecs = eigsh(H, k=2, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]
        E0_scan[k] = evals[0]
        gap_scan[k] = evals[1] - evals[0]
        sx_scan[k] = expectation(evecs[:, 0].astype(np.complex128), Sx0).real

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(h_scan, sx_scan, 'o-', color='crimson', lw=1.5, ms=5)
    for h0 in h_list:
        axes[0].axvline(h0, color='gray', ls=':', alpha=0.7)
    axes[0].axvline(1.0, color='black', ls='--', alpha=0.5, label='critical h_c = 1 (thermo. limit)')
    axes[0].set_xlabel("h")
    axes[0].set_ylabel(r"$\langle\sigma_1^x\rangle_{\mathrm{GS}}$")
    axes[0].set_title(f"Transverse magnetization vs h  (L={L})")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].plot(h_scan, gap_scan, 's-', color='darkgreen', lw=1.5, ms=5)
    for h0 in h_list:
        axes[1].axvline(h0, color='gray', ls=':', alpha=0.7)
    axes[1].axvline(1.0, color='black', ls='--', alpha=0.5)
    axes[1].set_xlabel("h")
    axes[1].set_ylabel(r"$E_1 - E_0$")
    axes[1].set_title(f"Energy gap vs h  (L={L})")
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    out2 = os.path.join(ASSET_DIR, "b1_scan_h.png")
    fig.savefig(out2, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"[Saved] {out2}")

    return results


if __name__ == "__main__":
    run_b1()
