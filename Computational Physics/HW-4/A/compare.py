"""
非简谐振子 vs 简谐振子 多维度对比分析
简谐振子解析解: E_n = n + 1/2, phi_n(x) = Hermite 函数
非简谐振子: V(x) = x^2/2 + x^4/10, 由 FEM 数值求解
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import hermite
from math import factorial, sqrt, pi
import matplotlib.pyplot as plt
import os

# ============================================================
# 1. FEM 求解非简谐振子（复用 solve.py 的逻辑）
# ============================================================
L = 10.0
N = 1000
num_states = 5

x_nodes = np.linspace(-L, L, N + 1)
h = x_nodes[1] - x_nodes[0]
N_inner = N - 1

def V_anharmonic(x):
    return 0.5 * x**2 + 0.1 * x**4

def V_harmonic(x):
    return 0.5 * x**2

def build_fem_matrices(V_func, x_nodes, N, N_inner):
    """组装 FEM 刚度矩阵和质量矩阵"""
    gauss_pts = np.array([-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)])
    gauss_wts = np.array([5.0/9.0, 8.0/9.0, 5.0/9.0])

    K = np.zeros((N_inner, N_inner))
    M = np.zeros((N_inner, N_inner))

    for e in range(N):
        x_left = x_nodes[e]
        x_right = x_nodes[e + 1]
        he = x_right - x_left

        K_kinetic = np.array([[1.0, -1.0], [-1.0, 1.0]]) / (2.0 * he)
        K_potential = np.zeros((2, 2))
        M_local = np.zeros((2, 2))

        for gp in range(3):
            xi = gauss_pts[gp]
            x_phys = 0.5 * (x_left + x_right) + 0.5 * he * xi
            w = gauss_wts[gp] * he / 2.0
            N1 = (x_right - x_phys) / he
            N2 = (x_phys - x_left) / he
            N_vec = np.array([N1, N2])
            K_potential += w * V_func(x_phys) * np.outer(N_vec, N_vec)
            M_local += w * np.outer(N_vec, N_vec)

        K_elem = K_kinetic + K_potential
        for i_local in range(2):
            i_global = e + i_local - 1
            if i_global < 0 or i_global >= N_inner:
                continue
            for j_local in range(2):
                j_global = e + j_local - 1
                if j_global < 0 or j_global >= N_inner:
                    continue
                K[i_global, j_global] += K_elem[i_local, j_local]
                M[i_global, j_global] += M_local[i_local, j_local]
    return K, M

def solve_fem(V_func, x_nodes, N, N_inner, num_states):
    """FEM 求解并返回能量和归一化波函数"""
    K, M = build_fem_matrices(V_func, x_nodes, N, N_inner)
    eigenvalues, eigenvectors = eigh(K, M, subset_by_index=[0, num_states - 1])
    phi_all = np.zeros((N + 1, num_states))
    phi_all[1:N, :] = eigenvectors
    for k in range(num_states):
        norm = np.trapz(phi_all[:, k]**2, x_nodes)
        phi_all[:, k] /= np.sqrt(norm)
    return eigenvalues, phi_all

# 求解两个系统
E_ah, phi_ah = solve_fem(V_anharmonic, x_nodes, N, N_inner, num_states)
E_ho, phi_ho = solve_fem(V_harmonic, x_nodes, N, N_inner, num_states)

# 简谐振子解析解
E_ho_exact = np.array([n + 0.5 for n in range(num_states)])

def harmonic_wavefunction(n, x):
    """简谐振子第 n 个本征态的解析波函数"""
    Hn = hermite(n)
    coeff = 1.0 / sqrt(2**n * factorial(n)) * (1.0 / pi)**0.25
    return coeff * np.exp(-x**2 / 2.0) * Hn(x)

phi_ho_exact = np.zeros((len(x_nodes), num_states))
for n in range(num_states):
    phi_ho_exact[:, n] = harmonic_wavefunction(n, x_nodes)

# ============================================================
# 输出路径
# ============================================================
asset_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "asset")
os.makedirs(asset_dir, exist_ok=True)

plt.rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

# ============================================================
# 对比 1: 能级对比
# ============================================================
print("=" * 70)
print("对比 1: 能级对比")
print("=" * 70)
print(f"{'态':>6} | {'简谐(解析)':>12} | {'简谐(FEM)':>12} | {'非简谐(FEM)':>12} | {'偏移量':>10} | {'偏移率':>10}")
print("-" * 70)
for k in range(num_states):
    label = f"n={k}"
    delta_E = E_ah[k] - E_ho_exact[k]
    ratio = delta_E / E_ho_exact[k] * 100
    print(f"{label:>6} | {E_ho_exact[k]:>12.6f} | {E_ho[k]:>12.6f} | {E_ah[k]:>12.6f} | {delta_E:>+10.6f} | {ratio:>+9.3f}%")
print("=" * 70)

# 能级间距对比
print("\n" + "=" * 70)
print("对比 2: 能级间距 Delta_E = E_{n+1} - E_n")
print("=" * 70)
print(f"{'间距':>10} | {'简谐(解析)':>12} | {'非简谐(FEM)':>12} | {'差异':>10}")
print("-" * 70)
for k in range(num_states - 1):
    dE_ho = E_ho_exact[k+1] - E_ho_exact[k]
    dE_ah = E_ah[k+1] - E_ah[k]
    diff = dE_ah - dE_ho
    print(f"E_{k+1}-E_{k}  | {dE_ho:>12.6f} | {dE_ah:>12.6f} | {diff:>+10.6f}")
print("=" * 70)
print("简谐振子能级等间距 (Delta_E = 1)；非简谐振子能级间距递增。")

# ============================================================
# 对比 3: 波函数统计量（<x^2>, <x^4>, 不确定度）
# ============================================================
print("\n" + "=" * 70)
print("对比 3: 期望值与不确定度")
print("=" * 70)

# 计算期望值
def compute_expectation(phi, x, op_func):
    """计算 <phi|op|phi>"""
    integrand = phi * op_func(x) * phi
    return np.trapz(integrand, x)

print(f"{'态':>6} | {'系统':>6} | {'<x^2>':>10} | {'<x^4>':>10} | {'Delta_x':>10} | {'<V>':>10} | {'<T>':>10}")
print("-" * 80)
for k in range(num_states):
    for label, phi_arr, E_arr in [("简谐", phi_ho_exact, E_ho_exact),
                                   ("非简谐", phi_ah, E_ah)]:
        phi_k = phi_arr[:, k]
        x2 = compute_expectation(phi_k, x_nodes, lambda x: x**2)
        x4 = compute_expectation(phi_k, x_nodes, lambda x: x**4)
        dx = np.sqrt(x2)  # <x>=0 by symmetry
        if label == "简谐":
            V_exp = compute_expectation(phi_k, x_nodes, V_harmonic)
        else:
            V_exp = compute_expectation(phi_k, x_nodes, V_anharmonic)
        T_exp = E_arr[k] - V_exp  # E = <T> + <V>
        print(f"n={k:>3} | {label:>6} | {x2:>10.6f} | {x4:>10.6f} | {dx:>10.6f} | {V_exp:>10.6f} | {T_exp:>10.6f}")
    print("-" * 80)

# ============================================================
# 对比 4: FEM 数值精度验证（简谐振子 FEM vs 解析解）
# ============================================================
print("\n" + "=" * 70)
print("对比 4: FEM 数值精度验证 (简谐振子 FEM vs 解析解)")
print("=" * 70)
print(f"{'态':>6} | {'E_exact':>12} | {'E_FEM':>12} | {'绝对误差':>12} | {'相对误差':>12}")
print("-" * 60)
for k in range(num_states):
    err_abs = abs(E_ho[k] - E_ho_exact[k])
    err_rel = err_abs / E_ho_exact[k]
    print(f"n={k:>3} | {E_ho_exact[k]:>12.6f} | {E_ho[k]:>12.6f} | {err_abs:>12.2e} | {err_rel:>12.2e}")
print("=" * 70)

# ============================================================
# 绘图
# ============================================================

# --- 图1: 能级对比柱状图 ---
fig, ax = plt.subplots(figsize=(10, 6))
x_pos = np.arange(num_states)
width = 0.35
bars1 = ax.bar(x_pos - width/2, E_ho_exact, width, label='简谐振子 (解析)', color='#4ECDC4', edgecolor='black')
bars2 = ax.bar(x_pos + width/2, E_ah, width, label='非简谐振子 (FEM)', color='#FF6B6B', edgecolor='black')
for i, (v1, v2) in enumerate(zip(E_ho_exact, E_ah)):
    ax.text(i - width/2, v1 + 0.1, f'{v1:.4f}', ha='center', fontsize=9)
    ax.text(i + width/2, v2 + 0.1, f'{v2:.4f}', ha='center', fontsize=9)
ax.set_xlabel('量子数 n', fontsize=14)
ax.set_ylabel('能量 E', fontsize=14)
ax.set_title('能级对比: 简谐振子 vs 非简谐振子', fontsize=16)
ax.set_xticks(x_pos)
ax.set_xticklabels([f'n={i}' for i in range(num_states)])
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_energy_comparison.png"))
print("\n已保存: asset/A_energy_comparison.png")

# --- 图2: 能级间距对比 ---
fig, ax = plt.subplots(figsize=(8, 5))
spacing_ho = np.diff(E_ho_exact)
spacing_ah = np.diff(E_ah)
x_pos = np.arange(len(spacing_ho))
ax.plot(x_pos, spacing_ho, 'o-', markersize=10, linewidth=2, label='简谐振子', color='#4ECDC4')
ax.plot(x_pos, spacing_ah, 's-', markersize=10, linewidth=2, label='非简谐振子', color='#FF6B6B')
ax.set_xlabel('间距编号', fontsize=14)
ax.set_ylabel(r'$\Delta E = E_{n+1} - E_n$', fontsize=14)
ax.set_title('能级间距对比', fontsize=16)
ax.set_xticks(x_pos)
ax.set_xticklabels([f'$E_{i+1}-E_{i}$' for i in range(len(spacing_ho))])
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_energy_spacing.png"))
print("已保存: asset/A_energy_spacing.png")

# --- 图3: 波函数逐态对比 ---
fig, axes = plt.subplots(num_states, 1, figsize=(10, 14), sharex=True)
for k in range(num_states):
    ax = axes[k]
    # 符号对齐
    sign = np.sign(np.dot(phi_ah[:, k], phi_ho_exact[:, k]))
    phi_ah_k = sign * phi_ah[:, k]
    ax.plot(x_nodes, phi_ho_exact[:, k], '--', linewidth=2, color='#4ECDC4',
            label='简谐 (解析)')
    ax.plot(x_nodes, phi_ah_k, '-', linewidth=1.5, color='#FF6B6B',
            label='非简谐 (FEM)')
    ax.set_xlim(-5, 5)
    ax.set_ylabel(r'$\phi_%d(x)$' % k, fontsize=12)
    ax.set_title(f'n={k}:  简谐 E={E_ho_exact[k]:.4f}  vs  非简谐 E={E_ah[k]:.4f}', fontsize=12)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
axes[-1].set_xlabel('x', fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_wavefunction_comparison.png"))
print("已保存: asset/A_wavefunction_comparison.png")

# --- 图4: 概率密度对比 ---
fig, axes = plt.subplots(num_states, 1, figsize=(10, 14), sharex=True)
for k in range(num_states):
    ax = axes[k]
    ax.plot(x_nodes, phi_ho_exact[:, k]**2, '--', linewidth=2, color='#4ECDC4',
            label='简谐')
    ax.plot(x_nodes, phi_ah[:, k]**2, '-', linewidth=1.5, color='#FF6B6B',
            label='非简谐')
    ax.fill_between(x_nodes, phi_ho_exact[:, k]**2, phi_ah[:, k]**2,
                     alpha=0.2, color='gray', label='差异区域')
    ax.set_xlim(-5, 5)
    ax.set_ylabel(r'$|\phi_%d|^2$' % k, fontsize=12)
    ax.set_title(f'n={k}: 概率密度对比', fontsize=12)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)
axes[-1].set_xlabel('x', fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_density_comparison.png"))
print("已保存: asset/A_density_comparison.png")

# --- 图5: 势能对比 ---
fig, ax = plt.subplots(figsize=(9, 6))
x_fine = np.linspace(-4, 4, 500)
ax.plot(x_fine, V_harmonic(x_fine), '--', linewidth=2, color='#4ECDC4',
        label=r'$V_{HO} = \frac{1}{2}x^2$')
ax.plot(x_fine, V_anharmonic(x_fine), '-', linewidth=2, color='#FF6B6B',
        label=r'$V_{AH} = \frac{1}{2}x^2 + \frac{1}{10}x^4$')
ax.fill_between(x_fine, V_harmonic(x_fine), V_anharmonic(x_fine),
                alpha=0.2, color='gray')
# 画能级
for k in range(num_states):
    ax.axhline(y=E_ho_exact[k], color='#4ECDC4', linewidth=0.8,
               linestyle=':', alpha=0.7)
    ax.axhline(y=E_ah[k], color='#FF6B6B', linewidth=0.8,
               linestyle=':', alpha=0.7)
ax.set_xlabel('x', fontsize=14)
ax.set_ylabel('V(x) / E', fontsize=14)
ax.set_title('势能与能级对比', fontsize=16)
ax.set_ylim(-0.5, 8)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_potential_comparison.png"))
print("已保存: asset/A_potential_comparison.png")

# --- 图6: 能级偏移量随量子数的变化 ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
delta_E = E_ah - E_ho_exact
ratio_E = delta_E / E_ho_exact * 100

ax1.bar(range(num_states), delta_E, color='#FF6B6B', edgecolor='black')
ax1.set_xlabel('量子数 n', fontsize=14)
ax1.set_ylabel(r'$\Delta E = E_{AH} - E_{HO}$', fontsize=14)
ax1.set_title('能级偏移量', fontsize=15)
ax1.set_xticks(range(num_states))
for i, v in enumerate(delta_E):
    ax1.text(i, v + 0.02, f'{v:.4f}', ha='center', fontsize=10)
ax1.grid(True, alpha=0.3, axis='y')

ax2.bar(range(num_states), ratio_E, color='#4ECDC4', edgecolor='black')
ax2.set_xlabel('量子数 n', fontsize=14)
ax2.set_ylabel('偏移率 (%)', fontsize=14)
ax2.set_title('能级偏移率', fontsize=15)
ax2.set_xticks(range(num_states))
for i, v in enumerate(ratio_E):
    ax2.text(i, v + 0.2, f'{v:.2f}%', ha='center', fontsize=10)
ax2.grid(True, alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_energy_shift.png"))
print("已保存: asset/A_energy_shift.png")

plt.close('all')
print("\n对比分析完成。")
