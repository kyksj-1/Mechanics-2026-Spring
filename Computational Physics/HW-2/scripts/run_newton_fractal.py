# -*- coding: utf-8 -*-
"""
B. 牛顿迭代法 — Newton 分形可视化脚本

功能:
    1. B.3: 以 (0,0) 为中心，半宽 1，分辨率 0.002 的 Newton 分形
    2. B.4a: 以 (-0.8, 0.0) 为中心，半宽 0.25，分辨率 0.0005 的放大
    3. B.4b: 以 (-0.56, 0.18) 为中心，半宽 0.1，分辨率 0.0002 的放大
    4. 额外: 迭代次数热力图、收敛速度分析、单点迭代轨迹可视化

使用方式:
    conda activate research_env
    python scripts/run_newton_fractal.py

作者: kyksj-1
"""

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")  # 非交互式后端
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from src.newton import (
    compute_newton_fractal,
    newton_iterate_general,
    ROOTS_Z3,
)
from config.config import OUTPUT_DIR

# ============================================================
# 配置 matplotlib 中文显示
# ============================================================
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# ============================================================
# 三个根对应的颜色方案
# ============================================================
# 根 0 (z=1): 红色, 根 1 (exp(2πi/3)): 绿色, 根 2 (exp(4πi/3)): 蓝色
# 未收敛: 黑色
ROOT_COLORS = ["#E74C3C", "#2ECC71", "#3498DB"]  # 红、绿、蓝
UNCONVERGED_COLOR = "#1A1A2E"                      # 深黑色


def make_root_colormap() -> ListedColormap:
    """
    构建离散颜色映射: -1 -> 黑色, 0 -> 红, 1 -> 绿, 2 -> 蓝。

    由于 imshow 需要连续 colormap，我们将 root_index 映射为 0,1,2,3
    (0=未收敛, 1=根0, 2=根1, 3=根2)。
    """
    return ListedColormap([UNCONVERGED_COLOR] + ROOT_COLORS)


def plot_fractal(
    root_index: np.ndarray,
    iter_count: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    title: str,
    save_name: str,
    show_iter: bool = True,
) -> None:
    """
    绘制 Newton 分形图: 根归属 + 迭代次数着色。

    Parameters
    ----------
    root_index : np.ndarray
        根编号数组 (ny, nx)，值为 -1, 0, 1, 2。
    iter_count : np.ndarray
        迭代次数数组 (ny, nx)。
    x_axis, y_axis : np.ndarray
        坐标轴。
    title : str
        图标题。
    save_name : str
        保存文件名（不含路径）。
    show_iter : bool
        是否同时显示迭代次数热力图。
    """
    # 将 root_index 从 {-1, 0, 1, 2} 映射到 {0, 1, 2, 3}
    display_data = root_index + 1  # -1 -> 0, 0 -> 1, 1 -> 2, 2 -> 3
    cmap = make_root_colormap()

    extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]

    if show_iter:
        fig, axes = plt.subplots(1, 2, figsize=(18, 8))
        fig.suptitle(title, fontsize=16, fontweight="bold")

        # --- 左图: 根归属 ---
        ax1 = axes[0]
        im1 = ax1.imshow(
            display_data,
            extent=extent,
            origin="lower",
            cmap=cmap,
            vmin=0, vmax=3,
            interpolation="nearest",
            aspect="equal",
        )
        ax1.set_xlabel("Re(z)", fontsize=12)
        ax1.set_ylabel("Im(z)", fontsize=12)
        ax1.set_title("收敛到的根")

        # 在图上标注三个根的位置
        for idx, root in enumerate(ROOTS_Z3):
            ax1.plot(
                root.real, root.imag,
                marker="*", markersize=12,
                color="white", markeredgecolor="black", markeredgewidth=1.0,
            )
            ax1.annotate(
                f"r{idx}",
                xy=(root.real, root.imag),
                xytext=(root.real + 0.05, root.imag + 0.05),
                fontsize=10, color="white",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.6),
            )

        # 添加颜色图例
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=ROOT_COLORS[0], label=f"r0 = 1"),
            Patch(facecolor=ROOT_COLORS[1], label=f"r1 = e^(2pi*i/3)"),
            Patch(facecolor=ROOT_COLORS[2], label=f"r2 = e^(4pi*i/3)"),
            Patch(facecolor=UNCONVERGED_COLOR, edgecolor="white", label="未收敛"),
        ]
        ax1.legend(handles=legend_elements, loc="upper right", fontsize=9)

        # --- 右图: 迭代次数热力图 ---
        ax2 = axes[1]
        im2 = ax2.imshow(
            iter_count,
            extent=extent,
            origin="lower",
            cmap="inferno",
            interpolation="nearest",
            aspect="equal",
        )
        ax2.set_xlabel("Re(z)", fontsize=12)
        ax2.set_ylabel("Im(z)", fontsize=12)
        ax2.set_title("收敛迭代次数")
        fig.colorbar(im2, ax=ax2, label="迭代次数", shrink=0.8)

    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(10, 10))
        fig.suptitle(title, fontsize=16, fontweight="bold")

        im1 = ax1.imshow(
            display_data,
            extent=extent,
            origin="lower",
            cmap=cmap,
            vmin=0, vmax=3,
            interpolation="nearest",
            aspect="equal",
        )
        ax1.set_xlabel("Re(z)", fontsize=12)
        ax1.set_ylabel("Im(z)", fontsize=12)

    plt.tight_layout()
    output_path = OUTPUT_DIR / save_name
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"  [保存] {output_path}")
    plt.close(fig)


def plot_shaded_fractal(
    root_index: np.ndarray,
    iter_count: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    title: str,
    save_name: str,
) -> None:
    """
    绘制根归属 + 迭代次数混合着色的高质量分形图。

    颜色由根编号决定色相，迭代次数决定亮度，
    迭代越多越暗，展现分形边界的精细结构。
    """
    # 构造 RGB 图像
    ny, nx = root_index.shape
    img = np.zeros((ny, nx, 3), dtype=np.float64)

    # 将颜色字符串解析为 RGB
    rgb_colors = []
    for c in ROOT_COLORS:
        r = int(c[1:3], 16) / 255.0
        g = int(c[3:5], 16) / 255.0
        b = int(c[5:7], 16) / 255.0
        rgb_colors.append((r, g, b))

    # 迭代次数归一化 (用于亮度调制)
    max_iter_val = float(iter_count.max()) if iter_count.max() > 0 else 1.0
    # 使用对数映射让亮度变化更明显
    brightness = 1.0 - 0.7 * (np.log1p(iter_count) / np.log1p(max_iter_val))

    # 为每个根着色
    for idx in range(3):
        mask = root_index == idx
        for ch in range(3):
            img[mask, ch] = rgb_colors[idx][ch] * brightness[mask]

    # 未收敛点设为黑色 (默认已为 0)

    extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]

    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.imshow(
        img,
        extent=extent,
        origin="lower",
        interpolation="nearest",
        aspect="equal",
    )
    ax.set_xlabel("Re(z)", fontsize=12)
    ax.set_ylabel("Im(z)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    plt.tight_layout()
    output_path = OUTPUT_DIR / save_name
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"  [保存] {output_path}")
    plt.close(fig)


def plot_iteration_trajectory(save_name: str) -> None:
    """
    绘制几个典型初始点的牛顿迭代轨迹，展示收敛过程。
    """
    # 选取几个有代表性的初始点
    starting_points = [
        0.5 + 0.5j,    # 应该收敛到 r0 或 r1
        -0.5 + 0.1j,   # 靠近 r1/r2 边界
        0.1 + 0.8j,    # 远离实轴
        -0.8 + 0.0j,   # 分形边界附近
    ]

    f = lambda z: z**3 - 1
    f_prime = lambda z: 3 * z**2

    fig, ax = plt.subplots(figsize=(10, 10))

    # 先画一个淡化的分形背景
    root_idx_bg, _, x_bg, y_bg = compute_newton_fractal(
        center=(0.0, 0.0), half_width=1.2, resolution=0.005, max_iter=50
    )
    display_bg = root_idx_bg + 1
    cmap_bg = make_root_colormap()
    ax.imshow(
        display_bg, extent=[x_bg[0], x_bg[-1], y_bg[0], y_bg[-1]],
        origin="lower", cmap=cmap_bg, vmin=0, vmax=3,
        interpolation="nearest", aspect="equal", alpha=0.3,
    )

    # 画三个根
    for idx, root in enumerate(ROOTS_Z3):
        ax.plot(root.real, root.imag, "*", markersize=15,
                color=ROOT_COLORS[idx], markeredgecolor="black", markeredgewidth=1.0)

    # 对每个初始点画迭代轨迹
    traj_colors = ["#F39C12", "#8E44AD", "#1ABC9C", "#E91E63"]
    for i, z0 in enumerate(starting_points):
        trajectory = [z0]
        z = z0
        for _ in range(30):
            fz = f(z)
            fpz = f_prime(z)
            if abs(fpz) < 1e-15:
                break
            z_new = z - fz / fpz
            trajectory.append(z_new)
            if abs(z_new - z) < 1e-10:
                break
            z = z_new

        traj = np.array(trajectory)
        ax.plot(traj.real, traj.imag, "o-", color=traj_colors[i],
                markersize=4, linewidth=1.5, alpha=0.9,
                label=f"z0 = {z0:.1f}")
        # 标记起点
        ax.plot(z0.real, z0.imag, "D", color=traj_colors[i],
                markersize=8, markeredgecolor="black")

    ax.set_xlabel("Re(z)", fontsize=12)
    ax.set_ylabel("Im(z)", fontsize=12)
    ax.set_title("牛顿迭代轨迹示例", fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    output_path = OUTPUT_DIR / save_name
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"  [保存] {output_path}")
    plt.close(fig)


def main() -> None:
    """牛顿分形可视化主流程"""

    print("=" * 60)
    print("B. 牛顿迭代法 — Newton 分形可视化")
    print("=" * 60)

    # --- 打印三个根 ---
    print("\nf(z) = z^3 - 1 的三个根 (三次单位根):")
    for i, root in enumerate(ROOTS_Z3):
        print(f"  r{i} = {root.real:+.6f} {root.imag:+.6f}i"
              f"  (|r{i}| = {abs(root):.6f})")

    # ==========================================================
    # B.3: 全局视图
    # ==========================================================
    print("\n--- B.3: 全局 Newton 分形 ---")
    print("  区域: [-1, 1] x [-1, 1]，分辨率 0.002")

    t0 = time.perf_counter()
    root_idx_3, iter_cnt_3, x3, y3 = compute_newton_fractal(
        center=(0.0, 0.0), half_width=1.0, resolution=0.002, max_iter=100
    )
    t1 = time.perf_counter()
    print(f"  计算耗时: {t1 - t0:.2f} s")
    print(f"  网格大小: {len(x3)} x {len(y3)} = {len(x3) * len(y3)} 个点")
    print(f"  未收敛点数: {np.sum(root_idx_3 == -1)}")

    # 基本分形图 (根归属 + 迭代次数)
    plot_fractal(root_idx_3, iter_cnt_3, x3, y3,
                 "B.3 Newton 分形 — 全局视图 (中心=(0,0), 半宽=1, 分辨率=0.002)",
                 "b3_newton_fractal_global.png")

    # 高质量着色版本 (亮度 = 迭代次数)
    plot_shaded_fractal(root_idx_3, iter_cnt_3, x3, y3,
                        "B.3 Newton 分形 (迭代次数着色)",
                        "b3_newton_fractal_shaded.png")

    # ==========================================================
    # B.4a: 放大到 (-0.8, 0.0) 附近
    # ==========================================================
    print("\n--- B.4a: 放大 — 中心 (-0.8, 0.0) ---")
    print("  区域: [-1.05, -0.55] x [-0.25, 0.25]，分辨率 0.0005")

    t0 = time.perf_counter()
    root_idx_4a, iter_cnt_4a, x4a, y4a = compute_newton_fractal(
        center=(-0.8, 0.0), half_width=0.25, resolution=0.0005, max_iter=200
    )
    t1 = time.perf_counter()
    print(f"  计算耗时: {t1 - t0:.2f} s")
    print(f"  网格大小: {len(x4a)} x {len(y4a)} = {len(x4a) * len(y4a)} 个点")

    plot_fractal(root_idx_4a, iter_cnt_4a, x4a, y4a,
                 "B.4a Newton 分形 — 放大 (中心=(-0.8, 0.0), 半宽=0.25, 分辨率=0.0005)",
                 "b4a_newton_fractal_zoom1.png")

    plot_shaded_fractal(root_idx_4a, iter_cnt_4a, x4a, y4a,
                        "B.4a Newton 分形 (迭代次数着色)",
                        "b4a_newton_fractal_zoom1_shaded.png")

    # ==========================================================
    # B.4b: 放大到 (-0.56, 0.18) 附近
    # ==========================================================
    print("\n--- B.4b: 放大 — 中心 (-0.56, 0.18) ---")
    print("  区域: [-0.66, -0.46] x [0.08, 0.28]，分辨率 0.0002")

    t0 = time.perf_counter()
    root_idx_4b, iter_cnt_4b, x4b, y4b = compute_newton_fractal(
        center=(-0.56, 0.18), half_width=0.1, resolution=0.0002, max_iter=300
    )
    t1 = time.perf_counter()
    print(f"  计算耗时: {t1 - t0:.2f} s")
    print(f"  网格大小: {len(x4b)} x {len(y4b)} = {len(x4b) * len(y4b)} 个点")

    plot_fractal(root_idx_4b, iter_cnt_4b, x4b, y4b,
                 "B.4b Newton 分形 — 放大 (中心=(-0.56, 0.18), 半宽=0.1, 分辨率=0.0002)",
                 "b4b_newton_fractal_zoom2.png")

    plot_shaded_fractal(root_idx_4b, iter_cnt_4b, x4b, y4b,
                        "B.4b Newton 分形 (迭代次数着色)",
                        "b4b_newton_fractal_zoom2_shaded.png")

    # ==========================================================
    # 额外可视化: 迭代轨迹
    # ==========================================================
    print("\n--- 额外: 牛顿迭代轨迹 ---")
    plot_iteration_trajectory("b_extra_iteration_trajectory.png")

    print("\n全部可视化完成!")


if __name__ == "__main__":
    main()
