"""
梳子晶格随机行走: 计算 ⟨|x|⟩ vs t, 与 1D 随机行走对比
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
from comb_walk import simulate_comb_walk, simulate_1d_walk

print("=" * 60)
print("梳子晶格上的随机行走")
print("=" * 60)

# ============================================================
# 运行参数
# ============================================================
num_walkers = 50000
max_steps = 10000
print(f"行走者数: {num_walkers}, 最大步数: {max_steps}")

# 梳子行走
print("运行梳子行走...")
x_comb, y_comb = simulate_comb_walk(num_walkers, max_steps, seed=42)
abs_x_comb = np.abs(x_comb)
mean_abs_x_comb = abs_x_comb.mean(axis=0)   # shape: (max_steps+1,)

# 1D 行走
print("运行 1D 行走...")
x_1d = simulate_1d_walk(num_walkers, max_steps, seed=123)
mean_abs_x_1d = np.abs(x_1d).mean(axis=0)

step_array = np.arange(max_steps + 1)

# ============================================================
# 输出结果
# ============================================================
print(f"\n关键数据点 (步数 → 梳子⟨|x|⟩ vs 1D⟨|x|⟩):")
checkpoints = [100, 500, 1000, 2000, 5000, 10000]
for cp in checkpoints:
    if cp <= max_steps:
        print(f"  t={cp:5d}  梳子⟨|x|⟩={mean_abs_x_comb[cp]:.4f}   1D⟨|x|⟩={mean_abs_x_1d[cp]:.4f}")

# ============================================================
# 拟合: log-log 斜率
# ============================================================
# 只在 t>=10 范围拟合以避免 t=0 问题
fit_start = 10
log_t = np.log(step_array[fit_start:])
log_comb = np.log(mean_abs_x_comb[fit_start:])
log_1d = np.log(mean_abs_x_1d[fit_start:])

slope_comb, intercept_comb = np.polyfit(log_t, log_comb, 1)
slope_1d, intercept_1d = np.polyfit(log_t, log_1d, 1)

print(f"\nlog-log 拟合斜率 (t≥{fit_start}):")
print(f"  梳子: ⟨|x|⟩ ∝ t^{slope_comb:.4f}  (预期: 1/4 = 0.25)")
print(f"  1D:   ⟨|x|⟩ ∝ t^{slope_1d:.4f}  (预期: 1/2 = 0.5)")
print(f"  梳子幂律常数: exp(intercept) = {np.exp(intercept_comb):.4f}")
print(f"  1D 幂律常数:   exp(intercept) = {np.exp(intercept_1d):.4f}")

# ============================================================
# 绘图1: ⟨|x|⟩ vs t (log-log)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 线性坐标
axes[0].plot(step_array, mean_abs_x_comb, "r-", label="梳子随机行走", linewidth=1)
axes[0].plot(step_array, mean_abs_x_1d, "b-", label="1D 随机行走", linewidth=1)
axes[0].set_xlabel("t (步数)")
axes[0].set_ylabel(r"$\langle |x| \rangle$")
axes[0].set_title(r"$\langle |x| \rangle$ vs t (线性坐标)")
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# log-log 坐标
axes[1].loglog(step_array[1:], mean_abs_x_comb[1:], "r-", label=f"梳子 (斜率 ≈ {slope_comb:.3f})", linewidth=1)
axes[1].loglog(step_array[1:], mean_abs_x_1d[1:], "b-", label=f"1D (斜率 ≈ {slope_1d:.3f})", linewidth=1)
# 参考线
t_ref = np.array([1, max_steps])
axes[1].loglog(t_ref, intercept_comb * t_ref**0.25, "r--", alpha=0.4, label=r"$t^{1/4}$")
axes[1].loglog(t_ref, intercept_1d * t_ref**0.5, "b--", alpha=0.4, label=r"$t^{1/2}$")
axes[1].set_xlabel("t (步数)")
axes[1].set_ylabel(r"$\langle |x| \rangle$")
axes[1].set_title(r"$\langle |x| \rangle$ vs t (log-log)")
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
asset_dir = os.path.join(os.path.dirname(__file__), "..", "asset")
os.makedirs(asset_dir, exist_ok=True)
fig.savefig(os.path.join(asset_dir, "fig_b_combwalk.png"), dpi=150, bbox_inches="tight")
plt.close()

# ============================================================
# 绘图2: 行走者空间分布快照
# ============================================================
# 取前 5000 个行走者在几个时刻的位置
fig2, axes2 = plt.subplots(2, 3, figsize=(14, 9))
snapshots = [10, 100, 1000, 3000, 7000, 10000]
sample_n = min(5000, num_walkers)

for idx, snap in enumerate(snapshots):
    ax = axes2[idx // 3, idx % 3]
    ax.scatter(x_comb[:sample_n, snap], y_comb[:sample_n, snap],
               s=0.5, c="red", alpha=0.5, rasterized=True)
    ax.set_xlim(-60, 60)
    ax.set_ylim(-60, 60)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"t = {snap}")
    ax.set_aspect("equal")
    ax.axhline(y=0, color="blue", linewidth=0.5, alpha=0.3)
    ax.axvline(x=0, color="blue", linewidth=0.5, alpha=0.3)
    ax.grid(True, alpha=0.2)

plt.suptitle("梳子随机行走: 不同时刻的粒子分布", fontsize=14)
plt.tight_layout()
fig2.savefig(os.path.join(asset_dir, "fig_b_snapshots.png"), dpi=150, bbox_inches="tight")
plt.close()

# ============================================================
# 绘图3: y分布 (验证: 大多数时间在主干附近)
# ============================================================
fig3, axes3 = plt.subplots(1, 2, figsize=(12, 5))

# y 方均根
rms_y = np.sqrt(np.mean(y_comb**2, axis=0))
rms_x_comb = np.sqrt(np.mean(x_comb**2, axis=0))

axes3[0].plot(step_array, rms_x_comb, "r-", label=r"$\sqrt{\langle x^2 \rangle}$ 梳子")
axes3[0].plot(step_array, rms_y, "g-", label=r"$\sqrt{\langle y^2 \rangle}$ 梳子")
axes3[0].set_xlabel("t")
axes3[0].set_ylabel("RMS 位移")
axes3[0].set_title("x 和 y 方向的方均根位移")
axes3[0].legend()
axes3[0].grid(True, alpha=0.3)

# 在主干上的概率比例 vs t
on_spine_prob = np.mean(y_comb == 0, axis=0)
axes3[1].plot(step_array, on_spine_prob, "b-", linewidth=1)
axes3[1].set_xlabel("t")
axes3[1].set_ylabel("P(y=0)")
axes3[1].set_title("行走者在主干上的概率 vs t")
axes3[1].grid(True, alpha=0.3)

plt.tight_layout()
fig3.savefig(os.path.join(asset_dir, "fig_b_rms.png"), dpi=150, bbox_inches="tight")
plt.close()

print("\n图片已保存: fig_b_combwalk.png, fig_b_snapshots.png, fig_b_rms.png")
print("完成!")
