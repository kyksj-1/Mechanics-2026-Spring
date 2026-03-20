# -*- coding: utf-8 -*-
"""
A.1 DFT 测试脚本

功能:
    1. 生成随机 complex128 数组
    2. 用自实现的 DFT 计算变换结果
    3. 与 numpy.fft.fft 的结果进行对比
    4. 可视化对比结果

使用方式:
    conda activate research_env
    python scripts/run_dft_test.py

作者: kyksj-1
"""

import sys
from pathlib import Path

# 将项目根目录加入 sys.path，以便导入 src 和 config
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")  # 非交互式后端，避免 GUI 阻塞
import matplotlib.pyplot as plt
from src.dft import dft, idft
from config.config import OUTPUT_DIR

# ============================================================
# 配置 matplotlib 中文显示
# ============================================================
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def main() -> None:
    """DFT 测试主流程"""

    print("=" * 60)
    print("A.1 DFT 测试: 自实现 DFT vs numpy.fft.fft")
    print("=" * 60)

    # ----------------------------------------------------------
    # 1. 生成测试数据: 随机 complex128 数组
    # ----------------------------------------------------------
    rng = np.random.default_rng(seed=42)  # 可复现的随机数生成器
    test_sizes = [8, 16, 32, 64, 128]     # 多个测试长度

    for N in test_sizes:
        # 生成随机复数数组: 实部和虚部均为标准正态分布
        x: np.ndarray = (
            rng.standard_normal(N) + 1j * rng.standard_normal(N)
        ).astype(np.complex128)

        # ----------------------------------------------------------
        # 2. 分别计算 DFT
        # ----------------------------------------------------------
        X_my: np.ndarray = dft(x)              # 自实现的 DFT
        X_np: np.ndarray = np.fft.fft(x)       # numpy 标准库

        # ----------------------------------------------------------
        # 3. 计算误差
        # ----------------------------------------------------------
        # 最大绝对误差
        max_abs_err: float = np.max(np.abs(X_my - X_np))
        # 相对误差 (L2 范数)
        rel_err: float = np.linalg.norm(X_my - X_np) / np.linalg.norm(X_np)

        print(f"\n--- N = {N} ---")
        print(f"  最大绝对误差: {max_abs_err:.2e}")
        print(f"  相对误差 (L2): {rel_err:.2e}")
        print(f"  通过: {'YES' if max_abs_err < 1e-10 else 'NO'}")

    # ----------------------------------------------------------
    # 4. 可视化: 取 N=64 的案例详细展示
    # ----------------------------------------------------------
    print("\n[可视化] 生成 N=64 的 DFT 对比图...")

    N_plot = 64
    x_plot = (
        rng.standard_normal(N_plot) + 1j * rng.standard_normal(N_plot)
    ).astype(np.complex128)

    X_my_plot = dft(x_plot)
    X_np_plot = np.fft.fft(x_plot)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"DFT 对比验证 (N={N_plot})", fontsize=16, fontweight="bold")

    k = np.arange(N_plot)

    # --- 子图1: 幅值谱对比 ---
    axes[0, 0].stem(k, np.abs(X_my_plot), linefmt="b-", markerfmt="bo",
                     basefmt="b-", label="自实现 DFT")
    axes[0, 0].stem(k, np.abs(X_np_plot), linefmt="r--", markerfmt="rx",
                     basefmt="r-", label="numpy.fft.fft")
    axes[0, 0].set_title("幅值谱 |X[k]|")
    axes[0, 0].set_xlabel("频率索引 k")
    axes[0, 0].set_ylabel("|X[k]|")
    axes[0, 0].legend()

    # --- 子图2: 相位谱对比 ---
    axes[0, 1].stem(k, np.angle(X_my_plot), linefmt="b-", markerfmt="bo",
                     basefmt="b-", label="自实现 DFT")
    axes[0, 1].stem(k, np.angle(X_np_plot), linefmt="r--", markerfmt="rx",
                     basefmt="r-", label="numpy.fft.fft")
    axes[0, 1].set_title("相位谱 arg(X[k])")
    axes[0, 1].set_xlabel("频率索引 k")
    axes[0, 1].set_ylabel("相位 (rad)")
    axes[0, 1].legend()

    # --- 子图3: 绝对误差 ---
    abs_err = np.abs(X_my_plot - X_np_plot)
    axes[1, 0].stem(k, abs_err, linefmt="g-", markerfmt="go", basefmt="g-")
    axes[1, 0].set_title("逐点绝对误差 |X_my - X_np|")
    axes[1, 0].set_xlabel("频率索引 k")
    axes[1, 0].set_ylabel("绝对误差")
    axes[1, 0].set_yscale("log")

    # --- 子图4: IDFT 还原验证 ---
    x_recovered = idft(X_my_plot)
    recovery_err = np.abs(x_plot - x_recovered)
    axes[1, 1].stem(np.arange(N_plot), recovery_err, linefmt="m-",
                     markerfmt="mo", basefmt="m-")
    axes[1, 1].set_title("IDFT 还原误差 |x - IDFT(DFT(x))|")
    axes[1, 1].set_xlabel("样本索引 n")
    axes[1, 1].set_ylabel("绝对误差")
    axes[1, 1].set_yscale("log")

    plt.tight_layout()

    # 保存图片
    output_path = OUTPUT_DIR / "a1_dft_comparison.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"[保存] 图片已保存至: {output_path}")
    plt.show()


if __name__ == "__main__":
    main()
