"""
A-Q5: 打靶法 (Shooting Method) 求解非简谐振子
与有限元法形成 "局部方法 vs 全局方法" 的对比

原理:
  将定态薛定谔方程视为初值问题 (IVP):
    phi''(x) = 2(V(x) - E) phi(x)
  从 x = -L 开始，给定 phi(-L) ~ 0, phi'(-L) 的初始条件，
  向右积分到 x = +L。调整能量 E 使 phi(+L) = 0 (边界条件)。
  利用二分法或 Brent 方法在 E 上搜索零点。
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import os

# ============================================================
# 参数
# ============================================================
L = 10.0           # 计算域 [-L, L]
num_states = 5     # 求解前 5 个态

def V(x):
    return 0.5 * x**2 + 0.1 * x**4

def schrodinger_ode(x, y, E):
    """
    定态薛定谔方程改写为一阶ODE组:
    y[0] = phi(x)
    y[1] = phi'(x)
    dy[0]/dx = y[1]
    dy[1]/dx = 2*(V(x) - E) * y[0]
    """
    phi, dphi = y
    d2phi = 2.0 * (V(x) - E) * phi
    return [dphi, d2phi]

def shoot(E, parity):
    """
    从 x=0 向右积分到 x=L，利用奇偶对称性:
    偶态: phi(0)=1, phi'(0)=0
    奇态: phi(0)=0, phi'(0)=1
    返回 phi(L) 作为匹配条件
    """
    if parity == 'even':
        y0 = [1.0, 0.0]
    else:
        y0 = [0.0, 1.0]

    sol = solve_ivp(schrodinger_ode, [0, L], y0, args=(E,),
                    method='RK45', rtol=1e-12, atol=1e-14,
                    dense_output=True)
    return sol.y[0, -1]

# ============================================================
# 搜索本征能量
# ============================================================
print("=" * 60)
print("打靶法 (Shooting Method) 求解非简谐振子")
print("V(x) = x^2/2 + x^4/10")
print("=" * 60)

# 根据量子数的奇偶性确定对称性
# n=0,2,4: 偶态; n=1,3: 奇态
eigenvalues_shoot = []
E_search_ranges = [
    (0.3, 0.8, 'even'),    # n=0
    (1.3, 2.2, 'odd'),     # n=1
    (2.5, 3.5, 'even'),    # n=2
    (4.0, 5.2, 'odd'),     # n=3
    (5.5, 7.0, 'even'),    # n=4
]

for n, (E_low, E_high, parity) in enumerate(E_search_ranges):
    # 用 Brent 方法找零点
    try:
        E_n = brentq(lambda E: shoot(E, parity), E_low, E_high,
                     xtol=1e-12, rtol=1e-14)
        eigenvalues_shoot.append(E_n)
        label = "基态" if n == 0 else f"第{n}激发态"
        print(f"E_{n} ({label}) = {E_n:.6f}  [parity: {parity}]")
    except ValueError as e:
        print(f"E_{n}: 在区间 [{E_low}, {E_high}] 未找到零点: {e}")

eigenvalues_shoot = np.array(eigenvalues_shoot)

# ============================================================
# 计算波函数
# ============================================================
x_plot = np.linspace(0, L, 2000)
wavefunctions = []

for n, E_n in enumerate(eigenvalues_shoot):
    parity = 'even' if n % 2 == 0 else 'odd'
    if parity == 'even':
        y0 = [1.0, 0.0]
    else:
        y0 = [0.0, 1.0]

    sol = solve_ivp(schrodinger_ode, [0, L], y0, args=(E_n,),
                    method='RK45', rtol=1e-12, atol=1e-14,
                    t_eval=x_plot)

    # 构造完整波函数（利用对称性）
    x_full = np.concatenate([-x_plot[::-1], x_plot[1:]])
    if parity == 'even':
        phi_full = np.concatenate([sol.y[0, ::-1], sol.y[0, 1:]])
    else:
        phi_full = np.concatenate([-sol.y[0, ::-1], sol.y[0, 1:]])

    # 归一化
    norm = np.trapz(phi_full**2, x_full)
    phi_full /= np.sqrt(norm)
    wavefunctions.append((x_full, phi_full))

# ============================================================
# 与 FEM 结果交叉验证
# ============================================================
# FEM 结果（从 solve.py 的输出）
E_fem = np.array([0.559164, 1.769608, 3.138942, 4.629573, 6.221558])

print("\n" + "=" * 60)
print("交叉验证: 打靶法 vs 有限元法")
print("=" * 60)
print(f"{'态':>6} | {'打靶法':>12} | {'有限元法':>12} | {'绝对差异':>12}")
print("-" * 60)
for k in range(num_states):
    diff = abs(eigenvalues_shoot[k] - E_fem[k])
    print(f"n={k:>3} | {eigenvalues_shoot[k]:>12.6f} | {E_fem[k]:>12.6f} | {diff:>12.2e}")
print("=" * 60)

# ============================================================
# 方法对比总结
# ============================================================
print("\n" + "=" * 60)
print("方法对比总结")
print("=" * 60)
print("""
+------------------+---------------------------+---------------------------+
|     对比维度      |      有限元法 (FEM)        |      打靶法 (Shooting)      |
+------------------+---------------------------+---------------------------+
| 方法类型         | 全局方法                   | 局部方法                    |
| 核心思想         | 变分原理 + 基函数展开       | 初值问题 + 边界匹配         |
| 求解方式         | 广义特征值问题 Kφ=EMφ      | ODE积分 + 根搜索            |
| 一次性输出       | 所有本征态                 | 单个本征态                  |
| 精度控制         | 网格加密(h-refinement)     | ODE积分精度 + 根搜索精度    |
| 对称性利用       | 不需要（自动满足）          | 需要手动分类奇偶态          |
| 高激发态         | 矩阵尺寸增大但仍可行       | 搜索区间需要预估            |
| 计算复杂度       | O(N^2) 存储 + 特征值分解   | 每个态 O(N) 但需多次搜索    |
| 推广到高维       | 自然推广                   | 困难                       |
| 实现复杂度       | 中等（需组装矩阵）         | 较低（只需 ODE solver）     |
+------------------+---------------------------+---------------------------+
""")

# ============================================================
# 绘图
# ============================================================
asset_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "asset")
os.makedirs(asset_dir, exist_ok=True)

plt.rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

# --- 打靶法波函数 ---
fig, axes = plt.subplots(num_states, 1, figsize=(10, 14), sharex=True)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
for k in range(num_states):
    ax = axes[k]
    x_full, phi_full = wavefunctions[k]
    label_str = "基态" if k == 0 else f"第{k}激发态"
    ax.plot(x_full, phi_full, color=colors[k], linewidth=1.5)
    ax.fill_between(x_full, phi_full, alpha=0.2, color=colors[k])
    ax.set_ylabel(r'$\phi_%d(x)$' % k, fontsize=12)
    ax.set_title(f'{label_str}: E_{k} = {eigenvalues_shoot[k]:.6f} (打靶法)', fontsize=13)
    ax.set_xlim(-6, 6)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
axes[-1].set_xlabel('x', fontsize=14)
fig.suptitle('打靶法求解结果', fontsize=16, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_shooting_eigenstates.png"))
print("\n已保存: asset/A_shooting_eigenstates.png")

# --- 两种方法能量对比图 ---
fig, ax = plt.subplots(figsize=(8, 5))
x_pos = np.arange(num_states)
width = 0.35
ax.bar(x_pos - width/2, E_fem, width, label='有限元法 (FEM)', color='#4ECDC4', edgecolor='black')
ax.bar(x_pos + width/2, eigenvalues_shoot, width, label='打靶法 (Shooting)', color='#FF6B6B', edgecolor='black')
ax.set_xlabel('量子数 n', fontsize=14)
ax.set_ylabel('能量 E', fontsize=14)
ax.set_title('有限元法 vs 打靶法: 能量对比', fontsize=16)
ax.set_xticks(x_pos)
ax.set_xticklabels([f'n={i}' for i in range(num_states)])
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(os.path.join(asset_dir, "A_fem_vs_shooting.png"))
print("已保存: asset/A_fem_vs_shooting.png")

plt.close('all')
print("\n打靶法求解完成。")
