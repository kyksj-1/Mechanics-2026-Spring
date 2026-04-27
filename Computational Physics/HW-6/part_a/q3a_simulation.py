"""
问题3(a): V = 1/2 x^2, T=1, 两个初态下的 Langevin 动力学系综平均

初态 A: q(0)=0, v(0)=1
初态 B: q(0)=4, v(0)=0
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import sys
import os

# ---- 中文字体配置 ----
_cjk = [f.name for f in fm.fontManager.ttflist if "SimHei" in f.name or "Microsoft YaHei" in f.name]
if _cjk:
    plt.rcParams["font.sans-serif"] = [_cjk[0]] + plt.rcParams.get("font.sans-serif", [])
    plt.rcParams["axes.unicode_minus"] = False

sys.path.insert(0, os.path.dirname(__file__))
from sde_solver import simulate_langevin, compute_ensemble_energies

# ============================================================
# 物理参数
# ============================================================
k_B, m, lam = 1.0, 1.0, 1.0
T = 1.0
D = np.sqrt(2.0 * lam * k_B * T)

# 数值参数
dt = 0.005
t_max = 10.0
num_steps = int(t_max / dt)
num_particles = 10000
t = np.linspace(0, t_max, num_steps + 1)

# 势能
def V(q):
    return 0.5 * q**2

def V_prime(q):
    return q

# ============================================================
# 初态 A: q(0)=0, v(0)=1
# ============================================================
print("=" * 50)
print("问题3(a): 初态 A -- q(0)=0, v(0)=1")
print(f"参数: k_B={k_B}, m={m}, λ={lam}, T={T}, 噪声强度 D²={D**2:.3f}")
print(f"数值: dt={dt}, 步数={num_steps}, 粒子数={num_particles}")
print(f"理论平衡值: ⟨E_k⟩={k_B*T/2}, ⟨V⟩={k_B*T/2}")
print("-" * 50)

np.random.seed(42)
q0 = np.zeros(num_particles)
v0 = np.ones(num_particles)
qA, vA = simulate_langevin(q0, v0, V_prime, m, lam, k_B, T, dt, num_steps, num_particles)
EkA, VA = compute_ensemble_energies(qA, vA, V, m)

# 输出弛豫后的平衡值 (取最后 20% 时间段的平均)
equil_idx_A = int(0.8 * len(t))
print(f"初态A 弛豫后 ⟨E_k⟩ = {EkA[equil_idx_A:].mean():.6f}")
print(f"初态A 弛豫后 ⟨V⟩   = {VA[equil_idx_A:].mean():.6f}")
print()

# ============================================================
# 初态 B: q(0)=4, v(0)=0
# ============================================================
print("=" * 50)
print("问题3(a): 初态 B -- q(0)=4, v(0)=0")
print("-" * 50)

np.random.seed(123)
q0 = np.full(num_particles, 4.0)
v0 = np.zeros(num_particles)
qB, vB = simulate_langevin(q0, v0, V_prime, m, lam, k_B, T, dt, num_steps, num_particles)
EkB, VB = compute_ensemble_energies(qB, vB, V, m)

equil_idx_B = int(0.8 * len(t))
print(f"初态B 弛豫后 ⟨E_k⟩ = {EkB[equil_idx_B:].mean():.6f}")
print(f"初态B 弛豫后 ⟨V⟩   = {VB[equil_idx_B:].mean():.6f}")
print()

# ============================================================
# 绘图
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 初态 A
axes[0].plot(t, EkA, "r-", label=r"$\langle E_k \rangle (t)$", linewidth=1.2)
axes[0].plot(t, VA, "b--", label=r"$\langle V \rangle (t)$", linewidth=1.2)
axes[0].axhline(y=k_B*T/2, color="gray", linestyle=":", alpha=0.6, label=f"理论 $k_B T/2 = {k_B*T/2}$")
axes[0].set_xlabel("t")
axes[0].set_ylabel("Energy")
axes[0].set_title(r"初态 A: $q(0)=0,\ v(0)=1$")
axes[0].legend(fontsize=9)
axes[0].set_xlim(0, 10)
axes[0].grid(True, alpha=0.3)

# 初态 B
axes[1].plot(t, EkB, "r-", label=r"$\langle E_k \rangle (t)$", linewidth=1.2)
axes[1].plot(t, VB, "b--", label=r"$\langle V \rangle (t)$", linewidth=1.2)
axes[1].axhline(y=k_B*T/2, color="gray", linestyle=":", alpha=0.6, label=f"理论 $k_B T/2 = {k_B*T/2}$")
axes[1].set_xlabel("t")
axes[1].set_ylabel("Energy")
axes[1].set_title(r"初态 B: $q(0)=4,\ v(0)=0$")
axes[1].legend(fontsize=9)
axes[1].set_xlim(0, 10)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()

asset_dir = os.path.join(os.path.dirname(__file__), "..", "asset")
os.makedirs(asset_dir, exist_ok=True)
save_path = os.path.join(asset_dir, "fig_a3a_energy.png")
fig.savefig(save_path, dpi=150, bbox_inches="tight")
print(f"图片已保存: {save_path}")
plt.close()

print("\n完成!")
