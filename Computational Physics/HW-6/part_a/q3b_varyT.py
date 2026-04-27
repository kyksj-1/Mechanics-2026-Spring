"""
问题3(b): V = 1/2 x^2, 改变温度, 计算充分弛豫后系统平均动能和平均势能

分析系综平均与时间平均的关系（各态历经性）
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

# 数值参数
dt = 0.005
t_max = 50.0          # 更长弛豫时间
num_steps = int(t_max / dt)
num_particles = 5000
t = np.linspace(0, t_max, num_steps + 1)

# 温度范围
T_list = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])

def V_harmonic(q):
    return 0.5 * q**2

def V_prime_harmonic(q):
    return q

print("=" * 60)
print("问题3(b): V = 1/2 x^2, 变温弛豫")
print(f"数值: dt={dt}, t_max={t_max}, 粒子数={num_particles}")
print(f"温度采样点: {T_list}")
print("-" * 60)

Ek_eq = np.zeros(len(T_list))
V_eq = np.zeros(len(T_list))
Ek_std = np.zeros(len(T_list))
V_std = np.zeros(len(T_list))

# 选择一个长轨迹的粒子做时间平均（验证各态历经）
time_avg_Ek = np.zeros(len(T_list))
time_avg_V = np.zeros(len(T_list))

equil_frac = 0.5  # 取后 50% 作为平衡态
equil_idx = int(equil_frac * len(t))

for i, T in enumerate(T_list):
    D = np.sqrt(2.0 * lam * k_B * T)

    # 系综模拟
    np.random.seed(42 + i)
    q0 = np.random.randn(num_particles) * 0.1
    v0 = np.random.randn(num_particles) * np.sqrt(k_B * T / m) * 0.1
    q, v = simulate_langevin(q0, v0, V_prime_harmonic, m, lam, k_B, T, dt, num_steps, num_particles)
    Ek, Vq = compute_ensemble_energies(q, v, V_harmonic, m)

    Ek_eq[i] = Ek[equil_idx:].mean()
    V_eq[i] = Vq[equil_idx:].mean()
    Ek_std[i] = Ek[equil_idx:].std()
    V_std[i] = Vq[equil_idx:].std()

    # 单粒子时间平均 (第0号粒子)
    Ek_time = 0.5 * m * v[0, equil_idx:]**2
    V_time = V_harmonic(q[0, equil_idx:])
    time_avg_Ek[i] = Ek_time.mean()
    time_avg_V[i] = V_time.mean()

    print(f"T={T:4.1f}  ⟨E_k⟩(∞)={Ek_eq[i]:.6f}±{Ek_std[i]:.6f}"
          f"  ⟨V⟩(∞)={V_eq[i]:.6f}±{V_std[i]:.6f}"
          f"  理论={k_B*T/2:.3f}"
          f"  时间平均E_k={time_avg_Ek[i]:.6f}  时间平均V={time_avg_V[i]:.6f}")

# ============================================================
# 绘图1: 系综平均 vs 温度
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].errorbar(T_list, Ek_eq, yerr=Ek_std, fmt="rs-", capsize=3, label=r"系综 $\langle E_k \rangle(\infty)$")
axes[0].plot(T_list, 0.5 * k_B * T_list, "k--", label=r"理论 $k_B T/2$")
axes[0].set_xlabel("T")
axes[0].set_ylabel(r"$\langle E_k \rangle(\infty)$")
axes[0].set_title("平衡态平均动能 vs 温度")
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].errorbar(T_list, V_eq, yerr=V_std, fmt="bs-", capsize=3, label=r"系综 $\langle V \rangle(\infty)$")
axes[1].plot(T_list, 0.5 * k_B * T_list, "k--", label=r"理论 $k_B T/2$")
axes[1].set_xlabel("T")
axes[1].set_ylabel(r"$\langle V \rangle(\infty)$")
axes[1].set_title("平衡态平均势能 vs 温度")
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
asset_dir = os.path.join(os.path.dirname(__file__), "..", "asset")
os.makedirs(asset_dir, exist_ok=True)
fig.savefig(os.path.join(asset_dir, "fig_a3b_harmonic.png"), dpi=150, bbox_inches="tight")
plt.close()

# ============================================================
# 绘图2: 系综平均 vs 时间平均（验证各态历经性）
# ============================================================
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

axes2[0].plot(T_list, Ek_eq, "r.-", label=r"系综平均 $\langle E_k \rangle_{\mathrm{ens}}$")
axes2[0].plot(T_list, time_avg_Ek, "b.--", label=r"时间平均 $\langle E_k \rangle_{\mathrm{time}}$")
axes2[0].plot(T_list, 0.5 * k_B * T_list, "k--", alpha=0.5, label=r"理论 $k_B T/2$")
axes2[0].set_xlabel("T")
axes2[0].set_ylabel("Energy")
axes2[0].set_title("动能: 系综平均 vs 时间平均")
axes2[0].legend()
axes2[0].grid(True, alpha=0.3)

axes2[1].plot(T_list, V_eq, "r.-", label=r"系综平均 $\langle V \rangle_{\mathrm{ens}}$")
axes2[1].plot(T_list, time_avg_V, "b.--", label=r"时间平均 $\langle V \rangle_{\mathrm{time}}$")
axes2[1].plot(T_list, 0.5 * k_B * T_list, "k--", alpha=0.5, label=r"理论 $k_B T/2$")
axes2[1].set_xlabel("T")
axes2[1].set_ylabel("Energy")
axes2[1].set_title("势能: 系综平均 vs 时间平均")
axes2[1].legend()
axes2[1].grid(True, alpha=0.3)

plt.tight_layout()
fig2.savefig(os.path.join(asset_dir, "fig_a3b_ergodicity.png"), dpi=150, bbox_inches="tight")
plt.close()

print("\n图片已保存: fig_a3b_harmonic.png, fig_a3b_ergodicity.png")
print("完成!")
