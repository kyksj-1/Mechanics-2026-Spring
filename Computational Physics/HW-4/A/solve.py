"""
非简谐振子定态薛定谔方程的有限元求解
势能: V(x) = x^2/2 + x^4/10
方法: Galerkin 有限元法 (线性基函数) -> 广义特征值问题 Kφ = EMφ
"""

import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import os

# ============================================================
# 参数设置
# ============================================================
L = 10.0          # 计算域 [-L, L]（足够大使波函数在边界处衰减到零）
N = 1000          # 单元数
num_states = 5    # 求解的本征态数目（基态 + 4个激发态）

# ============================================================
# 网格生成
# ============================================================
x_nodes = np.linspace(-L, L, N + 1)  # N+1 个节点
h = x_nodes[1] - x_nodes[0]          # 均匀网格步长
N_inner = N - 1                       # 内部节点数（去掉两端 Dirichlet 边界）

def V(x):
    """势能函数"""
    return 0.5 * x**2 + 0.1 * x**4

# ============================================================
# 组装刚度矩阵 K 和质量矩阵 M
# ============================================================
# 对于线性基函数, 每个单元 [x_e, x_{e+1}] 上:
#   刚度矩阵 K 的贡献:
#     K_local = (1/(2h)) * [[1, -1], [-1, 1]]   (动能项)
#     + 高斯积分计算 V(x) 的贡献               (势能项)
#   质量矩阵 M 的贡献:
#     M_local = (h/6) * [[2, 1], [1, 2]]

# 使用 2 点 Gauss 积分精确积分势能项（对 x^4 * 线性基 的乘积）
# 实际上 x^4 * N_i * N_j 最高为6次多项式，需要4点Gauss积分
# 使用 3 点 Gauss 积分（精确到5次多项式，对本问题足够精确）
gauss_pts_ref = np.array([-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)])
gauss_wts_ref = np.array([5.0/9.0, 8.0/9.0, 5.0/9.0])

K = np.zeros((N_inner, N_inner))
M = np.zeros((N_inner, N_inner))

for e in range(N):
    # 单元 e: 从节点 e 到节点 e+1
    x_left = x_nodes[e]
    x_right = x_nodes[e + 1]
    he = x_right - x_left  # 单元长度

    # 动能项: (1/2) * integral(dN_i/dx * dN_j/dx, x_left, x_right)
    # 对于线性基函数: dN1/dx = -1/he, dN2/dx = 1/he
    K_kinetic = np.array([[1.0, -1.0],
                          [-1.0, 1.0]]) / (2.0 * he)

    # 势能项和质量矩阵: 用 Gauss 积分
    K_potential = np.zeros((2, 2))
    M_local = np.zeros((2, 2))

    for gp in range(len(gauss_pts_ref)):
        # 参考坐标 -> 物理坐标
        xi = gauss_pts_ref[gp]
        x_phys = 0.5 * (x_left + x_right) + 0.5 * he * xi
        w = gauss_wts_ref[gp] * he / 2.0  # 积分权重 * Jacobian

        # 线性基函数在 Gauss 点的值
        N1 = (x_right - x_phys) / he
        N2 = (x_phys - x_left) / he
        N_vec = np.array([N1, N2])

        # 势能项贡献
        V_val = V(x_phys)
        K_potential += w * V_val * np.outer(N_vec, N_vec)

        # 质量矩阵贡献
        M_local += w * np.outer(N_vec, N_vec)

    K_elem = K_kinetic + K_potential

    # 组装到全局矩阵（仅内部节点，编号 0..N_inner-1 对应节点 1..N-1）
    for i_local in range(2):
        i_global = e + i_local - 1  # 全局内部节点编号
        if i_global < 0 or i_global >= N_inner:
            continue
        for j_local in range(2):
            j_global = e + j_local - 1
            if j_global < 0 or j_global >= N_inner:
                continue
            K[i_global, j_global] += K_elem[i_local, j_local]
            M[i_global, j_global] += M_local[i_local, j_local]

# ============================================================
# 求解广义特征值问题: K phi = E M phi
# ============================================================
eigenvalues, eigenvectors = eigh(K, M, subset_by_index=[0, num_states - 1])

# 内部节点的波函数值（边界为零）
phi_all = np.zeros((N + 1, num_states))
phi_all[1:N, :] = eigenvectors

# 归一化: integral |phi|^2 dx = 1 (使用梯形积分)
for k in range(num_states):
    norm = np.trapz(phi_all[:, k]**2, x_nodes)
    phi_all[:, k] /= np.sqrt(norm)

# ============================================================
# 输出结果
# ============================================================
print("=" * 60)
print("非简谐振子 V(x) = x^2/2 + x^4/10 的本征能量")
print("有限元法 (线性基函数, N=%d 个单元, 域 [-%.1f, %.1f])" % (N, L, L))
print("=" * 60)
for k in range(num_states):
    label = "基态" if k == 0 else "第%d激发态" % k
    print("E_%d (%s) = %.6f" % (k, label, eigenvalues[k]))
print("=" * 60)

# ============================================================
# 绘图
# ============================================================
asset_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "asset")
os.makedirs(asset_dir, exist_ok=True)

plt.rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

# --- 图1: 基态波函数绝对值 ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(x_nodes, np.abs(phi_all[:, 0]), 'b-', linewidth=2, label=r'$|\phi_0(x)|$')
ax.set_xlabel('x', fontsize=14)
ax.set_ylabel(r'$|\phi_0(x)|$', fontsize=14)
ax.set_title('基态波函数绝对值 (归一化)', fontsize=16)
ax.legend(fontsize=12)
ax.set_xlim(-5, 5)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_ground_state_wavefunction.png"))
print("已保存: asset/A_ground_state_wavefunction.png")

# --- 图2: 前5个本征态波函数 ---
fig, axes = plt.subplots(5, 1, figsize=(10, 14), sharex=True)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
for k in range(num_states):
    ax = axes[k]
    label = "基态" if k == 0 else "第%d激发态" % k
    ax.plot(x_nodes, phi_all[:, k], color=colors[k], linewidth=1.5)
    ax.fill_between(x_nodes, phi_all[:, k], alpha=0.2, color=colors[k])
    ax.set_ylabel(r'$\phi_%d(x)$' % k, fontsize=12)
    ax.set_title(r'%s: $E_%d = %.6f$' % (label, k, eigenvalues[k]), fontsize=13)
    ax.set_xlim(-6, 6)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
axes[-1].set_xlabel('x', fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_all_eigenstates.png"))
print("已保存: asset/A_all_eigenstates.png")

# --- 图3: 概率密度分布 ---
fig, ax = plt.subplots(figsize=(10, 6))
for k in range(num_states):
    label = r'$|\phi_%d|^2$, $E_%d=%.4f$' % (k, k, eigenvalues[k])
    ax.plot(x_nodes, phi_all[:, k]**2, color=colors[k], linewidth=1.5, label=label)
ax.set_xlabel('x', fontsize=14)
ax.set_ylabel(r'$|\phi(x)|^2$', fontsize=14)
ax.set_title('概率密度分布', fontsize=16)
ax.set_xlim(-6, 6)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_probability_density.png"))
print("已保存: asset/A_probability_density.png")

# --- 图4: 势能与能级示意图 ---
fig, ax = plt.subplots(figsize=(9, 6))
x_fine = np.linspace(-4, 4, 500)
ax.plot(x_fine, V(x_fine), 'k-', linewidth=2, label=r'$V(x) = \frac{1}{2}x^2 + \frac{1}{10}x^4$')
for k in range(num_states):
    # 经典转折点
    E_k = eigenvalues[k]
    # 画能级线
    # 找经典允许区域
    mask = V(x_fine) <= E_k
    if np.any(mask):
        x_class = x_fine[mask]
        x_left_tp = x_class[0]
        x_right_tp = x_class[-1]
        ax.hlines(E_k, x_left_tp, x_right_tp, colors=colors[k], linewidth=2)
        ax.text(x_right_tp + 0.2, E_k, r'$E_%d=%.4f$' % (k, E_k),
                fontsize=10, va='center', color=colors[k])
ax.set_xlabel('x', fontsize=14)
ax.set_ylabel('E', fontsize=14)
ax.set_title('势能曲线与能级图', fontsize=16)
ax.set_ylim(-0.5, eigenvalues[-1] + 2)
ax.legend(fontsize=12, loc='upper left')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_potential_energy_levels.png"))
print("已保存: asset/A_potential_energy_levels.png")

plt.close('all')
print("\nA 题求解完成。")
