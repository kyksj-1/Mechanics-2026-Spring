"""
Part B.2: Quench 动力学

t<0: h = 0.5 的 GS
t>=0: 哈密顿量瞬变为 h = 3.0
计算 <sigma_1^x(t)> 在 t in [0,20], 比较精确演化 (谱分解) 与 RK4.
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh

from tfim_core import build_tfim_hamiltonian, sigma_x_op, expectation


ASSET_DIR = os.path.join(os.path.dirname(__file__), "asset")


def evolve_exact(H, psi0, t_arr):
    """对 H 完整对角化, 用谱分解给出 |psi(t)> = sum c_n e^{-i E_n t} |n>.

    向量化: phase = exp(-i E_n t_k) * c_n  (Nt, N),  psi_t = phase @ evecs^T  (Nt, N)
    """
    H_dense = H.toarray()
    evals, evecs = np.linalg.eigh(H_dense)
    coeff = evecs.conj().T @ psi0  # (N,)
    # phase shape (Nt, N): outer
    phase = np.exp(-1j * np.outer(t_arr, evals)) * coeff[np.newaxis, :]
    psi_t = phase @ evecs.T  # (Nt, N)
    return psi_t


def evolve_rk4(H, psi0, t_arr, dt_max=0.01):
    """以 RK4 求解 i dpsi/dt = H psi, 即 dpsi/dt = -i H psi.

    在请求时间点 t_arr 之间用足够小的步长积分, 然后用 H 的稀疏矩乘.
    """
    def rhs(psi):
        return -1j * (H @ psi)

    psi_t = np.empty((len(t_arr), H.shape[0]), dtype=np.complex128)
    psi = psi0.astype(np.complex128).copy()
    psi_t[0] = psi
    for k in range(1, len(t_arr)):
        T_target = t_arr[k]
        t_now = t_arr[k - 1]
        # 选 dt: 不超过 dt_max, 且整除 (T_target - t_now)
        span = T_target - t_now
        nsub = max(1, int(np.ceil(span / dt_max)))
        dt = span / nsub
        for _ in range(nsub):
            k1 = rhs(psi)
            k2 = rhs(psi + 0.5 * dt * k1)
            k3 = rhs(psi + 0.5 * dt * k2)
            k4 = rhs(psi + dt * k3)
            psi = psi + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        psi_t[k] = psi
    return psi_t


def run_b2():
    L = 12
    J = 1.0
    h_i = 0.5  # 初始
    h_f = 3.0  # quench 后
    t_arr = np.linspace(0, 20, 401)  # 401 个时间点 (含端点)

    print("=" * 70)
    print(f"Part B.2: L={L}, J={J}, quench h: {h_i} -> {h_f},  t∈[0,20]")
    print("=" * 70)

    # 1) 求初态 (h_i 下 GS)
    H_i = build_tfim_hamiltonian(L, J, h_i, pbc=True)
    evals_i, evecs_i = eigsh(H_i, k=1, which='SA')
    psi0 = evecs_i[:, 0].astype(np.complex128)
    print(f"初态: h={h_i} 的基态, E_0 = {evals_i[0]:.8f}")

    # 2) 演化哈密顿量
    H_f = build_tfim_hamiltonian(L, J, h_f, pbc=True)
    Sx0 = sigma_x_op(L, site=0)

    # 精确演化
    t0 = time.time()
    psi_exact = evolve_exact(H_f, psi0, t_arr)
    sx_exact = np.array([expectation(psi_exact[k], Sx0).real for k in range(len(t_arr))])
    t_exact = time.time() - t0
    print(f"精确谱分解演化耗时: {t_exact:.3f} s")

    # RK4 演化
    t0 = time.time()
    psi_rk4 = evolve_rk4(H_f, psi0, t_arr, dt_max=0.005)
    sx_rk4 = np.array([expectation(psi_rk4[k], Sx0).real for k in range(len(t_arr))])
    t_rk = time.time() - t0
    print(f"RK4 (dt_max=0.005) 演化耗时: {t_rk:.3f} s")

    # 检查归一性
    norm_exact = np.linalg.norm(psi_exact, axis=1)
    norm_rk4 = np.linalg.norm(psi_rk4, axis=1)
    print(f"|psi|  精确:  min={norm_exact.min():.10f}, max={norm_exact.max():.10f}")
    print(f"|psi|  RK4 :  min={norm_rk4.min():.10f},  max={norm_rk4.max():.10f}")

    # 两种方法的差异
    abs_err = np.abs(sx_exact - sx_rk4)
    print(f"\n两种方法 <sigma_1^x>(t) 的差异: max = {abs_err.max():.3e}, mean = {abs_err.mean():.3e}")

    # 几个抽样输出
    print(f"\n抽样比较 (t, exact, RK4, |err|):")
    for tt in [0.0, 1.0, 5.0, 10.0, 15.0, 20.0]:
        idx = np.argmin(np.abs(t_arr - tt))
        print(f"  t={t_arr[idx]:5.2f}   exact = {sx_exact[idx]:+.10f}   RK4 = {sx_rk4[idx]:+.10f}   |Δ| = {abs(sx_exact[idx] - sx_rk4[idx]):.3e}")

    # 绘图
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                             gridspec_kw={'height_ratios': [3, 1]})
    axes[0].plot(t_arr, sx_exact, '-', color='royalblue', lw=2.0, label='Exact (spectral)')
    axes[0].plot(t_arr, sx_rk4, '--', color='crimson', lw=1.4, label='RK4 (dt=0.005)')
    axes[0].set_ylabel(r"$\langle \sigma_1^x(t) \rangle$", fontsize=13)
    axes[0].set_title(f"Quench dynamics: $h={h_i}\\rightarrow{h_f}$,  $L={L}$,  PBC")
    axes[0].legend(loc='upper right')
    axes[0].grid(alpha=0.3)

    axes[1].semilogy(t_arr, np.maximum(abs_err, 1e-16), color='black', lw=1.0)
    axes[1].set_xlabel("t")
    axes[1].set_ylabel("|exact - RK4|")
    axes[1].set_title("Pointwise discrepancy")
    axes[1].grid(alpha=0.3, which='both')

    fig.tight_layout()
    out = os.path.join(ASSET_DIR, "b2_quench_sigma_x.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\n[Saved] {out}")


if __name__ == "__main__":
    run_b2()
