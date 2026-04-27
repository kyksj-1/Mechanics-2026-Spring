"""
问题3(c): V = 1/2 x^4, 改变温度, 计算充分弛豫后系统平均动能和平均势能

与 3(b) 对比：四次方势不再满足简单的能均分定理
使用 scipy.integrate.quad 数值计算 Boltzmann 分布下的理论平均势能
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
from scipy.integrate import quad

# ============================================================
# 物理参数
# ============================================================
k_B, m, lam = 1.0, 1.0, 1.0

# 数值参数 — 缩减规模以保证运行效率
dt = 0.005
t_max = 30.0
num_steps = int(t_max / dt)
num_particles = 3000
t = np.linspace(0, t_max, num_steps + 1)

T_list = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])

def V_quartic(q):
    return 0.5 * q**4

def V_prime_quartic(q):
    return 2.0 * q**3

# ============================================================
# 理论值: V = q^4/2 的 Boltzmann 平均
# ============================================================
# 代入 y = q / (k_B T)^(1/4), 得到:
#   ⟨V⟩ = k_B T · ⟨y^4⟩_0 / 2
# 其中 ⟨y^4⟩_0 = ∫ y^4 exp(-y^4/2) dy / ∫ exp(-y^4/2) dy 与 T 无关
# 数值计算此无量纲比值:

def mean_y4_over_2():
    """计算无量纲量 ⟨y^4/2⟩_0, 其中权重 exp(-y^4/2)"""
    y_max = 10.0
    num, _ = quad(lambda y: y**4/2 * np.exp(-y**4/2), -y_max, y_max, limit=200)
    den, _ = quad(lambda y: np.exp(-y**4/2), -y_max, y_max, limit=200)
    return num / den

C_q4 = mean_y4_over_2()
print(f"无量纲常数 ⟨y^4/2⟩_0 = {C_q4:.6f}")
print(f"因此 ⟨V⟩ = {C_q4:.6f} · k_B T")

def theory_V_q4(T):
    return C_q4 * k_B * T

def theory_Ek(T):
    return k_B * T / 2.0

print("=" * 60)
print(f"问题3(c): V = (1/2) x^4, 变温弛豫")
print(f"数值: dt={dt}, t_max={t_max}, 步数={num_steps}, 粒子数={num_particles}")
print(f"温度采样点: {T_list}")
print(f"理论: ⟨E_k⟩ = k_B T/2,  ⟨V⟩ = {C_q4:.4f} · k_B T")
print("-" * 60)

Ek_eq = np.zeros(len(T_list))
V_eq = np.zeros(len(T_list))
Ek_std = np.zeros(len(T_list))
V_std = np.zeros(len(T_list))

equil_frac = 0.5
equil_idx = int(equil_frac * len(t))

for i, T in enumerate(T_list):
    np.random.seed(42 + i)
    # 小范围随机初始化
    q0 = np.random.randn(num_particles) * 0.1
    v0 = np.random.randn(num_particles) * np.sqrt(k_B * T / m) * 0.1
    q, v = simulate_langevin(q0, v0, V_prime_quartic, m, lam, k_B, T, dt, num_steps, num_particles)
    Ek, Vq = compute_ensemble_energies(q, v, V_quartic, m)

    Ek_eq[i] = Ek[equil_idx:].mean()
    V_eq[i] = Vq[equil_idx:].mean()
    Ek_std[i] = Ek[equil_idx:].std()
    V_std[i] = Vq[equil_idx:].std()

    print(f"T={T:4.1f}  ⟨E_k⟩={Ek_eq[i]:.6f}±{Ek_std[i]:.6f}"
          f"  ⟨V⟩={V_eq[i]:.6f}±{V_std[i]:.6f}"
          f"  Ek理论={theory_Ek(T):.4f}"
          f"  V理论={theory_V_q4(T):.6f}")

# ============================================================
# 绘图1: 动能和势能 vs 温度 (四次方势)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].errorbar(T_list, Ek_eq, yerr=Ek_std, fmt="rs-", capsize=3, label=r"$\langle E_k \rangle$ (系综)")
axes[0].plot(T_list, theory_Ek(T_list), "k--", label=r"$k_B T/2$ (能均分)")
axes[0].set_xlabel("T")
axes[0].set_ylabel(r"$\langle E_k \rangle(\infty)$")
axes[0].set_title("四次方势: 平衡态平均动能")
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].errorbar(T_list, V_eq, yerr=V_std, fmt="bs-", capsize=3, label=r"$\langle V \rangle$ (系综)")
axes[1].plot(T_list, theory_V_q4(T_list), "k--", label=f"理论 ${C_q4:.3f} k_B T$")
axes[1].plot(T_list, theory_Ek(T_list), "g:", alpha=0.5, label=r"$k_B T/2$ (参考)")
axes[1].set_xlabel("T")
axes[1].set_ylabel(r"$\langle V \rangle(\infty)$")
axes[1].set_title("四次方势: 平衡态平均势能")
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
asset_dir = os.path.join(os.path.dirname(__file__), "..", "asset")
os.makedirs(asset_dir, exist_ok=True)
fig.savefig(os.path.join(asset_dir, "fig_a3c_quartic.png"), dpi=150, bbox_inches="tight")
plt.close()

# ============================================================
# 绘图2: 二次 vs 四次 势能对比
# ============================================================
# 重新计算二次势数据以便在同时间刻度下对比
Ek_harm = np.zeros(len(T_list))
V_harm = np.zeros(len(T_list))

def V_harmonic(q):
    return 0.5 * q**2
def V_prime_harmonic(q):
    return q

for i, T in enumerate(T_list):
    np.random.seed(200 + i)
    q0 = np.random.randn(num_particles) * 0.1
    v0 = np.random.randn(num_particles) * np.sqrt(k_B * T / m) * 0.1
    q, v = simulate_langevin(q0, v0, V_prime_harmonic, m, lam, k_B, T, dt, num_steps, num_particles)
    Ek, Vq = compute_ensemble_energies(q, v, V_harmonic, m)
    Ek_harm[i] = Ek[equil_idx:].mean()
    V_harm[i] = Vq[equil_idx:].mean()

fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

axes2[0].plot(T_list, Ek_eq / T_list, "r.-", label="四次方势")
axes2[0].plot(T_list, Ek_harm / T_list, "b.-", label="二次方势")
axes2[0].axhline(y=k_B/2, color="gray", linestyle="--", alpha=0.5, label="$k_B/2$")
axes2[0].set_xlabel("T")
axes2[0].set_ylabel(r"$\langle E_k \rangle / T$")
axes2[0].set_title("动能/T vs T (二次 vs 四次)")
axes2[0].legend()
axes2[0].grid(True, alpha=0.3)

axes2[1].plot(T_list, V_eq / T_list, "r.-", label="四次方势")
axes2[1].plot(T_list, V_harm / T_list, "b.-", label="二次方势")
axes2[1].plot(T_list, theory_V_q4(T_list) / T_list, "r--", alpha=0.7, label="四次方理论")
axes2[1].axhline(y=k_B/2, color="gray", linestyle="--", alpha=0.5, label="$k_B/2$")
axes2[1].set_xlabel("T")
axes2[1].set_ylabel(r"$\langle V \rangle / T$")
axes2[1].set_title("势能/T vs T (二次 vs 四次)")
axes2[1].legend()
axes2[1].grid(True, alpha=0.3)

plt.tight_layout()
fig2.savefig(os.path.join(asset_dir, "fig_a3c_comparison.png"), dpi=150, bbox_inches="tight")
plt.close()

print("\n图片已保存: fig_a3c_quartic.png, fig_a3c_comparison.png")
print("完成!")
