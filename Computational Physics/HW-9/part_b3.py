"""
Part B.3: L=18, h=1, Néel 初态时间演化, 计算 Fidelity F(t) = |<psi(0)|psi(t)>|^2.

dim = 2^18 = 262144. 不能做完整对角化, 用 Krylov 子空间演化:
scipy.sparse.linalg.expm_multiply 在 start/stop/num 模式下一次性给出多个时间点.
为对照, 同时跑 RK4 验证一致性.
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import expm_multiply

from tfim_core import build_tfim_hamiltonian, neel_state_vector


ASSET_DIR = os.path.join(os.path.dirname(__file__), "asset")


def evolve_rk4_sparse(H, psi0, t_arr, dt_max=0.01):
    """RK4: dpsi/dt = -i H psi, H 稀疏."""
    psi_t = np.empty((len(t_arr), H.shape[0]), dtype=np.complex128)
    psi = psi0.astype(np.complex128).copy()
    psi_t[0] = psi
    for k in range(1, len(t_arr)):
        T_target = t_arr[k]
        t_now = t_arr[k - 1]
        span = T_target - t_now
        nsub = max(1, int(np.ceil(span / dt_max)))
        dt = span / nsub
        for _ in range(nsub):
            k1 = -1j * (H @ psi)
            k2 = -1j * (H @ (psi + 0.5 * dt * k1))
            k3 = -1j * (H @ (psi + 0.5 * dt * k2))
            k4 = -1j * (H @ (psi + dt * k3))
            psi = psi + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        psi_t[k] = psi
    return psi_t


def run_b3():
    L = 18
    J = 1.0
    h = 1.0
    N = 1 << L
    print("=" * 70)
    print(f"Part B.3: L={L}, J={J}, h={h}  (dim = {N})")
    print("初态: Néel |↑↓↑↓...> (题面 1-indexed 奇格点朝上)")
    print("=" * 70)

    t0 = time.time()
    H = build_tfim_hamiltonian(L, J, h, pbc=True)
    print(f"构造 H 耗时: {time.time() - t0:.3f} s,  nnz = {H.nnz}")

    psi0 = neel_state_vector(L)

    # 时间网格
    Nt = 401
    t_start, t_stop = 0.0, 10.0
    t_arr = np.linspace(t_start, t_stop, Nt)

    # --- 方法 1: Krylov (expm_multiply, 多时间点高效计算 exp(-iHt) psi0) ---
    print("\n[方法 1] Krylov (scipy expm_multiply, start/stop/num)...")
    t0 = time.time()
    # expm_multiply 计算 exp(s*A) @ B; 我们要 exp(-i*t*H) psi0, 即 A = -i*H, s=t.
    # 直接传入复矩阵 -1j*H, 它支持稀疏复阵.
    A = (-1j) * H
    # 注意: A 为反厄米, exp(s*A) 是酉矩阵, 数值上稳定.
    psi_t = expm_multiply(A, psi0, start=t_start, stop=t_stop, num=Nt, endpoint=True)
    print(f"  耗时: {time.time() - t0:.3f} s,  psi_t.shape = {psi_t.shape}")

    # Fidelity
    overlaps = psi_t @ psi0.conj()   # <psi(0)|psi(t)> = sum_i psi_t[k,i] * conj(psi0[i])
    # 注意: vdot(psi0, psi_t[k]) = sum conj(psi0)*psi_t[k]; 上面等价.
    fid = np.abs(overlaps) ** 2

    # 归一性检查
    norms = np.linalg.norm(psi_t, axis=1)
    print(f"  |psi(t)|: min={norms.min():.10f}, max={norms.max():.10f}")

    # --- 方法 2: RK4 (作为交叉验证, 取较小时间间隔) ---
    print("\n[方法 2] RK4 (dt_max=0.005, 用于验证)...")
    t0 = time.time()
    psi_rk4 = evolve_rk4_sparse(H, psi0, t_arr, dt_max=0.005)
    print(f"  耗时: {time.time() - t0:.3f} s")
    overlaps_rk = psi_rk4 @ psi0.conj()
    fid_rk = np.abs(overlaps_rk) ** 2
    norms_rk = np.linalg.norm(psi_rk4, axis=1)
    print(f"  |psi(t)| RK4: min={norms_rk.min():.10f}, max={norms_rk.max():.10f}")

    diff = np.abs(fid - fid_rk)
    print(f"\n两种方法 fidelity 差异: max={diff.max():.3e}, mean={diff.mean():.3e}")

    # 输出采样
    print(f"\n采样 (t, F_krylov, F_RK4):")
    for tt in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]:
        idx = np.argmin(np.abs(t_arr - tt))
        print(f"  t={t_arr[idx]:5.2f}   F_Krylov = {fid[idx]:.10f}   F_RK4 = {fid_rk[idx]:.10f}")

    # 画图
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                             gridspec_kw={'height_ratios': [3, 1]})
    axes[0].plot(t_arr, fid, '-', color='royalblue', lw=2.0, label='Krylov (expm_multiply)')
    axes[0].plot(t_arr, fid_rk, '--', color='crimson', lw=1.2, label='RK4 (dt=0.005)')
    axes[0].set_ylabel(r"$F(t) = |\langle\psi(0)|\psi(t)\rangle|^2$", fontsize=13)
    axes[0].set_title(f"Fidelity vs t: Néel quench under TFIM,  $L={L}$, $h={h}$, $J={J}$ (PBC)")
    axes[0].legend(loc='upper right')
    axes[0].grid(alpha=0.3)
    axes[0].set_yscale('log')
    axes[0].set_ylim(max(fid.min(), 1e-12) * 0.5, 1.5)

    axes[1].plot(t_arr, fid, '-', color='royalblue', lw=2.0)
    axes[1].set_xlabel("t")
    axes[1].set_ylabel("F(t) (linear)")
    axes[1].set_ylim(-0.02, 1.05)
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(ASSET_DIR, "b3_neel_fidelity.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\n[Saved] {out}")


if __name__ == "__main__":
    run_b3()
