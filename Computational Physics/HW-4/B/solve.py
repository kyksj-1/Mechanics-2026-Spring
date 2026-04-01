"""
B: 格点上的单粒子动力学
紧束缚模型 H = -sum(a_i^dag a_{i+1} + h.c.), 周期边界条件
L=200, i0=100, t in [0,50]

方法: 精确对角化
  H 在实空间为 L x L 循环三对角矩阵
  H 的本征值: E_k = -2 cos(2*pi*k/L), k=0,1,...,L-1
  时间演化: psi(t) = U exp(-iDt) U^dag psi(0)
"""

import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import os

# ============================================================
# 参数
# ============================================================
L = 200        # 格点数
i0 = 100       # 初始位置 (1-indexed, 程序中用 0-indexed: i0_idx = 99)
i0_idx = i0 - 1
t_max = 50.0   # 最大时间
t_snapshots = [1, 10, 20, 50]  # 绘图时刻

# ============================================================
# 构造哈密顿量 (周期边界条件的循环三对角矩阵)
# ============================================================
H = np.zeros((L, L))
for i in range(L):
    j_right = (i + 1) % L
    H[i, j_right] = -1.0
    H[j_right, i] = -1.0

# 对角化
eigenvalues, U = eigh(H)
# U: 列为本征向量, eigenvalues: 本征值
# 验证: E_k = -2*cos(2*pi*k/L)

# 初始波函数: psi_i(0) = delta(i - i0)
psi0 = np.zeros(L)
psi0[i0_idx] = 1.0

# 在本征基下展开: c_k = U^dag * psi0
c = U.T @ psi0

# ============================================================
# 时间演化函数
# ============================================================
def evolve(t):
    """精确时间演化 psi(t) = U * exp(-i*E*t) * c"""
    phase = np.exp(-1j * eigenvalues * t)
    psi_t = U @ (phase * c)
    return psi_t

# ============================================================
# Q1: 特定时刻的粒子密度分布
# ============================================================
print("=" * 60)
print("B: 格点上的单粒子动力学")
print(f"L={L}, i0={i0}, t in [0, {t_max}]")
print("=" * 60)

asset_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "asset")
os.makedirs(asset_dir, exist_ok=True)

plt.rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

sites = np.arange(1, L + 1)  # 格点编号 1 到 L

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes_flat = axes.flatten()

for idx, t in enumerate(t_snapshots):
    ax = axes_flat[idx]
    psi_t = evolve(t)
    rho = np.abs(psi_t)**2

    ax.plot(sites, rho, '-', linewidth=0.8, color='#1f77b4')
    ax.fill_between(sites, rho, alpha=0.3, color='#1f77b4')
    ax.axvline(x=i0, color='red', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_xlabel('格点编号 i', fontsize=12)
    ax.set_ylabel(r'$\rho_i(t) = |\psi_i(t)|^2$', fontsize=12)
    ax.set_title(f't = {t}', fontsize=14)
    ax.set_xlim(1, L)
    ax.grid(True, alpha=0.3)

    # 输出关键统计
    max_rho = np.max(rho)
    rho_i0 = rho[i0_idx]
    w = np.sqrt(np.sum((sites - i0)**2 * rho))
    print(f"t={t:>5}: rho_max={max_rho:.6f}, rho_i0={rho_i0:.6f}, w(t)={w:.4f}")

fig.suptitle(r'粒子密度分布 $\rho_i(t) = |\psi_i(t)|^2$', fontsize=16, y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "B_density_snapshots.png"), bbox_inches='tight')
print("\n已保存: asset/B_density_snapshots.png")

# ============================================================
# Q2: 波包宽度 w(t) 随时间的变化
# ============================================================
t_array = np.linspace(0, t_max, 2000)
w_array = np.zeros(len(t_array))

for idx, t in enumerate(t_array):
    psi_t = evolve(t)
    rho = np.abs(psi_t)**2
    w_array[idx] = np.sqrt(np.sum((sites - i0)**2 * rho))

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(t_array, w_array, '-', linewidth=1.5, color='#2ca02c')
ax.set_xlabel('时间 t', fontsize=14)
ax.set_ylabel(r'波包宽度 $w(t) = \sqrt{\sum_i (i-i_0)^2 \rho_i(t)}$', fontsize=13)
ax.set_title('波包宽度随时间的变化', fontsize=16)
ax.grid(True, alpha=0.3)

# 标注线性增长参考线
# 在 L -> inf 极限下, w(t) ~ sqrt(2) * t (群速度 v_max = 2)
t_ref = np.linspace(0.5, t_max, 100)
ax.plot(t_ref, np.sqrt(2) * t_ref, '--', color='gray', linewidth=1,
        label=r'$w = \sqrt{2}\,t$ (参考)')
ax.legend(fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "B_wavepacket_width.png"))
print("已保存: asset/B_wavepacket_width.png")

# ============================================================
# Q3: 起始点粒子密度 rho_{i0}(t) 随时间的变化
# ============================================================
rho_i0_array = np.zeros(len(t_array))

for idx, t in enumerate(t_array):
    psi_t = evolve(t)
    rho_i0_array[idx] = np.abs(psi_t[i0_idx])**2

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(t_array, rho_i0_array, '-', linewidth=1, color='#d62728')
ax.set_xlabel('时间 t', fontsize=14)
ax.set_ylabel(r'$\rho_{i_0}(t) = |\psi_{i_0}(t)|^2$', fontsize=14)
ax.set_title(r'起始点粒子密度 $\rho_{i_0}(t)$ 随时间的变化', fontsize=16)
ax.grid(True, alpha=0.3)

# 标注 1/t 衰减包络参考线
t_ref2 = np.linspace(1, t_max, 500)
# L->inf 时, rho_{i0}(t) = |J_0(2t)|^2 ~ 1/(pi*t) (渐近)
scale = 1.0 / (np.pi * t_ref2)
ax.plot(t_ref2, scale, '--', color='gray', linewidth=1,
        label=r'$\sim 1/(\pi t)$ 包络 (渐近)')
ax.legend(fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "B_origin_density.png"))
print("已保存: asset/B_origin_density.png")

# ============================================================
# 补充图: 时空密度演化图 (用于报告)
# ============================================================
t_fine = np.linspace(0, t_max, 500)
rho_spacetime = np.zeros((len(t_fine), L))

for idx, t in enumerate(t_fine):
    psi_t = evolve(t)
    rho_spacetime[idx, :] = np.abs(psi_t)**2

fig, ax = plt.subplots(figsize=(12, 6))
# 只显示中心区域
i_range = slice(30, 170)
im = ax.imshow(rho_spacetime[:, i_range].T, aspect='auto',
               extent=[0, t_max, 170, 31],
               cmap='hot', interpolation='bilinear')
ax.set_xlabel('时间 t', fontsize=14)
ax.set_ylabel('格点编号 i', fontsize=14)
ax.set_title('时空密度演化图', fontsize=16)
ax.axhline(y=i0, color='cyan', linewidth=0.8, linestyle='--', alpha=0.7)
plt.colorbar(im, ax=ax, label=r'$\rho_i(t)$')
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "B_spacetime_density.png"))
print("已保存: asset/B_spacetime_density.png")

# ============================================================
# 补充图: 色散关系
# ============================================================
fig, ax = plt.subplots(figsize=(8, 5))
k_values = np.linspace(-np.pi, np.pi, 500)
E_k = -2 * np.cos(k_values)
v_g = 2 * np.sin(k_values)  # 群速度 dE/dk

ax2 = ax.twinx()
ax.plot(k_values, E_k, '-', linewidth=2, color='#1f77b4', label=r'$E(k) = -2\cos k$')
ax2.plot(k_values, v_g, '--', linewidth=1.5, color='#ff7f0e', label=r'$v_g(k) = 2\sin k$')
ax.set_xlabel('波矢 k', fontsize=14)
ax.set_ylabel('能量 E(k)', fontsize=14, color='#1f77b4')
ax2.set_ylabel(r'群速度 $v_g(k)$', fontsize=14, color='#ff7f0e')
ax.set_title('紧束缚模型色散关系与群速度', fontsize=16)
ax.legend(fontsize=12, loc='upper left')
ax2.legend(fontsize=12, loc='upper right')
ax.grid(True, alpha=0.3)
ax.set_xlim(-np.pi, np.pi)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "B_dispersion.png"))
print("已保存: asset/B_dispersion.png")

plt.close('all')
print("\nB 题求解完成。")
