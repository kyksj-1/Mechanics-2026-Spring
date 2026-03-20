# -*- coding: utf-8 -*-
"""
A.2 FFT 基准测试脚本

功能:
    1. 对自实现的 Base-2 FFT 在 N=2^4 ~ 2^12 上进行正确性验证
    2. 精确计时对比自实现 FFT vs numpy.fft.fft
    3. 绘制计算时间 vs 数组大小的对比图，验证 O(N log N) 复杂度
    4. 与 O(N^2) DFT 进行复杂度对比

使用方式:
    conda activate research_env
    python scripts/run_fft_benchmark.py

作者: kyksj-1
"""

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")  # 非交互式后端，避免 GUI 阻塞
import matplotlib.pyplot as plt
import yaml
from src.fft import fft
from src.dft import dft
from config.config import OUTPUT_DIR

# ============================================================
# 配置 matplotlib 中文显示
# ============================================================
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def precise_timeit(func, *args, n_repeats: int = 5) -> float:
    """
    精确计时函数，返回多次运行的中位数耗时（秒）。

    使用 time.perf_counter_ns() 实现纳秒级精度。

    Parameters
    ----------
    func : callable
        要计时的函数。
    *args
        传给 func 的参数。
    n_repeats : int
        重复次数，取中位数以减少噪声。

    Returns
    -------
    float
        中位数耗时（秒）。
    """
    times = []
    for _ in range(n_repeats):
        # 使用 perf_counter_ns 获得纳秒级精度
        start_ns: int = time.perf_counter_ns()
        func(*args)
        end_ns: int = time.perf_counter_ns()
        times.append((end_ns - start_ns) * 1e-9)  # 转为秒
    return float(np.median(times))


def main() -> None:
    """FFT 基准测试主流程"""

    # ----------------------------------------------------------
    # 读取配置
    # ----------------------------------------------------------
    config_path = PROJECT_ROOT / "config" / "config.yaml"
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    bench_cfg = cfg["benchmark"]
    min_power: int = bench_cfg["min_power"]      # 最小指数
    max_power: int = bench_cfg["max_power"]      # 最大指数
    n_repeats: int = bench_cfg["n_repeats"]      # 重复次数
    seed: int = bench_cfg["seed"]                # 随机种子

    rng = np.random.default_rng(seed=seed)
    powers = list(range(min_power, max_power + 1))
    sizes = [2**p for p in powers]

    print("=" * 70)
    print("A.2 Base-2 FFT 基准测试")
    print("=" * 70)
    print(f"测试范围: 2^{min_power} ~ 2^{max_power} (N={sizes[0]} ~ {sizes[-1]})")
    print(f"每个大小重复 {n_repeats} 次，取中位数")
    print()

    # ----------------------------------------------------------
    # 1. 正确性验证 + 计时
    # ----------------------------------------------------------
    results = {
        "N": [],
        "my_fft_time": [],      # 自实现 FFT 耗时
        "numpy_fft_time": [],   # numpy FFT 耗时
        "my_dft_time": [],      # 自实现 DFT 耗时（仅小规模）
        "max_error": [],        # 最大绝对误差
    }

    # DFT 仅在 N <= 2^10 时测试（太大会很慢）
    dft_max_power = 10

    print(f"{'N':>8} | {'自实现FFT':>12} | {'numpy FFT':>12} | {'自实现DFT':>12} | {'最大误差':>12} | {'通过':>4}")
    print("-" * 78)

    for p in powers:
        N = 2**p
        # 生成随机复数测试数据
        x = (rng.standard_normal(N) + 1j * rng.standard_normal(N)).astype(np.complex128)

        # --- 正确性验证 ---
        X_my = fft(x)
        X_np = np.fft.fft(x)
        max_err = float(np.max(np.abs(X_my - X_np)))
        passed = max_err < 1e-8

        # --- 计时: 自实现 FFT ---
        t_my_fft = precise_timeit(fft, x, n_repeats=n_repeats)

        # --- 计时: numpy FFT ---
        t_np_fft = precise_timeit(np.fft.fft, x, n_repeats=n_repeats)

        # --- 计时: 自实现 DFT（仅小规模） ---
        if p <= dft_max_power:
            t_my_dft = precise_timeit(dft, x, n_repeats=n_repeats)
        else:
            t_my_dft = np.nan

        results["N"].append(N)
        results["my_fft_time"].append(t_my_fft)
        results["numpy_fft_time"].append(t_np_fft)
        results["my_dft_time"].append(t_my_dft)
        results["max_error"].append(max_err)

        # 格式化时间显示（自动选择单位）
        def fmt_time(t: float) -> str:
            if np.isnan(t):
                return "N/A"
            if t < 1e-3:
                return f"{t * 1e6:.1f} us"
            else:
                return f"{t * 1e3:.2f} ms"

        print(
            f"{N:>8} | {fmt_time(t_my_fft):>12} | {fmt_time(t_np_fft):>12} | "
            f"{fmt_time(t_my_dft):>12} | {max_err:>12.2e} | {'YES' if passed else 'NO':>4}"
        )

    # ----------------------------------------------------------
    # 2. 复杂度分析与可视化
    # ----------------------------------------------------------
    print("\n[可视化] 生成性能对比图...")

    N_arr = np.array(results["N"])
    t_my_fft = np.array(results["my_fft_time"])
    t_np_fft = np.array(results["numpy_fft_time"])
    t_my_dft = np.array(results["my_dft_time"])

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle("FFT 性能基准测试", fontsize=16, fontweight="bold")

    # --- 子图1: 绝对耗时对比 (双对数坐标) ---
    ax1 = axes[0]
    ax1.loglog(N_arr, t_my_fft, "o-", color="blue", linewidth=2,
               markersize=6, label="自实现 FFT")
    ax1.loglog(N_arr, t_np_fft, "s-", color="red", linewidth=2,
               markersize=6, label="numpy.fft.fft")

    # DFT 仅绘制有效数据点
    valid_dft = ~np.isnan(t_my_dft)
    ax1.loglog(N_arr[valid_dft], t_my_dft[valid_dft], "^--", color="green",
               linewidth=2, markersize=6, label="自实现 DFT (O(N^2))")

    # 叠加理论复杂度参考线
    # O(N log N) 参考线: 以自实现 FFT 的第一个点为锚点
    ref_nlogn = t_my_fft[0] * (N_arr * np.log2(N_arr)) / (N_arr[0] * np.log2(N_arr[0]))
    ax1.loglog(N_arr, ref_nlogn, ":", color="blue", alpha=0.4, label="O(N log N) 参考")

    # O(N^2) 参考线
    if valid_dft.any():
        first_valid = np.where(valid_dft)[0][0]
        ref_n2 = t_my_dft[first_valid] * (N_arr**2) / (N_arr[first_valid]**2)
        ax1.loglog(N_arr, ref_n2, ":", color="green", alpha=0.4, label="O(N^2) 参考")

    ax1.set_xlabel("数组大小 N", fontsize=12)
    ax1.set_ylabel("耗时 (秒)", fontsize=12)
    ax1.set_title("计算耗时 vs 数组大小")
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, which="both")

    # --- 子图2: 自实现 FFT vs numpy 的速度比 ---
    ax2 = axes[1]
    speed_ratio = t_my_fft / t_np_fft
    ax2.bar(
        [f"2^{p}" for p in powers],
        speed_ratio,
        color="steelblue",
        edgecolor="navy",
        alpha=0.8,
    )
    ax2.axhline(y=1.0, color="red", linestyle="--", linewidth=1.5, label="numpy 基准线")
    ax2.set_xlabel("数组大小", fontsize=12)
    ax2.set_ylabel("耗时比 (自实现 / numpy)", fontsize=12)
    ax2.set_title("相对 numpy 的速度比")
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis="y")

    # 在柱子上标注数值
    for i, ratio in enumerate(speed_ratio):
        ax2.text(i, ratio + 0.1, f"{ratio:.1f}x", ha="center", fontsize=9)

    plt.tight_layout()

    output_path = OUTPUT_DIR / "a2_fft_benchmark.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"[保存] 图片已保存至: {output_path}")
    plt.show()

    # ----------------------------------------------------------
    # 3. 复杂度理论分析输出
    # ----------------------------------------------------------
    print("\n" + "=" * 70)
    print("复杂度理论分析")
    print("=" * 70)
    print("""
自实现 Base-2 FFT 理论复杂度: O(N log N)

推导:
  - 递归关系: T(N) = 2·T(N/2) + O(N)
  - 每层递归做 N 次蝶形乘加运算
  - 共 log₂(N) 层递归
  - 总运算量: N·log₂(N) 次复数乘法和加法

对比 DFT 的 O(N²):
  - 当 N=4096 (2^12) 时，FFT 约需 49152 次运算
  - 同样大小的 DFT 需要 ~16,777,216 次运算
  - 加速比约 341 倍
""")


if __name__ == "__main__":
    main()
