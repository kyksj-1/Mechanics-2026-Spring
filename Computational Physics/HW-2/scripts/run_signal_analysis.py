# -*- coding: utf-8 -*-
"""
A.3 信号频域分析脚本

功能:
    1. 加载 waveform.dat 示波器数据
    2. 绘制时域波形
    3. 使用 numpy.fft 计算频谱
    4. 绘制频谱图，标注主要频率分量
    5. 可视化频谱的细节特征

使用方式:
    conda activate research_env
    python scripts/run_signal_analysis.py

作者: kyksj-1
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")  # 非交互式后端，避免 GUI 阻塞
import matplotlib.pyplot as plt
import yaml
from src.signal_analysis import (
    load_waveform,
    compute_frequency_spectrum,
    find_dominant_frequencies,
)
from config.config import OUTPUT_DIR, WAVEFORM_FILE

# ============================================================
# 配置 matplotlib 中文显示
# ============================================================
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def main() -> None:
    """信号频域分析主流程"""

    # ----------------------------------------------------------
    # 读取配置
    # ----------------------------------------------------------
    config_path = PROJECT_ROOT / "config" / "config.yaml"
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    sig_cfg = cfg["signal_analysis"]
    save_figures: bool = sig_cfg["save_figures"]
    figure_dpi: int = sig_cfg["figure_dpi"]

    print("=" * 60)
    print("A.3 信号频域分析")
    print("=" * 60)

    # ----------------------------------------------------------
    # 1. 加载数据
    # ----------------------------------------------------------
    print(f"\n[加载] 数据文件: {WAVEFORM_FILE}")
    time_arr, voltage_arr = load_waveform(str(WAVEFORM_FILE))

    N: int = len(time_arr)
    dt: float = time_arr[1] - time_arr[0]
    fs: float = 1.0 / dt
    T_total: float = time_arr[-1] - time_arr[0]

    print(f"  采样点数: {N}")
    print(f"  采样间隔: {dt:.6f} s")
    print(f"  采样率:   {fs:.1f} Hz")
    print(f"  总时长:   {T_total:.3f} s")
    print(f"  Nyquist 频率: {fs / 2:.1f} Hz")
    print(f"  电压范围: [{voltage_arr.min():.4f}, {voltage_arr.max():.4f}] V")

    # ----------------------------------------------------------
    # 2. 计算频谱
    # ----------------------------------------------------------
    print("\n[计算] FFT 频谱分析...")
    freqs, magnitude, spectrum = compute_frequency_spectrum(time_arr, voltage_arr)

    # ----------------------------------------------------------
    # 3. 找主要频率分量
    # ----------------------------------------------------------
    dom_freqs, dom_mags = find_dominant_frequencies(freqs, magnitude, threshold_ratio=0.05)

    print(f"\n[结果] 检测到 {len(dom_freqs)} 个主要频率分量:")
    print(f"{'频率 (Hz)':>12} | {'幅值':>10} | {'周期 (s)':>12}")
    print("-" * 42)
    for f_val, m_val in zip(dom_freqs, dom_mags):
        period = 1.0 / f_val if f_val > 0 else float("inf")
        print(f"{f_val:>12.2f} | {m_val:>10.4f} | {period:>12.6f}")

    # ----------------------------------------------------------
    # 4. 可视化: 综合分析图
    # ----------------------------------------------------------
    print("\n[可视化] 生成综合分析图...")

    fig = plt.figure(figsize=(18, 14))
    fig.suptitle("waveform.dat 信号频域分析", fontsize=18, fontweight="bold")

    # --- 子图1: 完整时域波形 ---
    ax1 = fig.add_subplot(3, 2, 1)
    ax1.plot(time_arr, voltage_arr, color="steelblue", linewidth=0.3, alpha=0.8)
    ax1.set_title("完整时域波形")
    ax1.set_xlabel("时间 (s)")
    ax1.set_ylabel("电压 (V)")
    ax1.grid(True, alpha=0.3)

    # --- 子图2: 时域波形局部放大 (前 0.5s) ---
    ax2 = fig.add_subplot(3, 2, 2)
    mask_zoom = time_arr <= 0.5
    ax2.plot(time_arr[mask_zoom], voltage_arr[mask_zoom],
             color="steelblue", linewidth=0.5)
    ax2.set_title("时域波形 (前 0.5s 局部放大)")
    ax2.set_xlabel("时间 (s)")
    ax2.set_ylabel("电压 (V)")
    ax2.grid(True, alpha=0.3)

    # --- 子图3: 完整频谱 (线性坐标) ---
    ax3 = fig.add_subplot(3, 2, 3)
    ax3.plot(freqs, magnitude, color="crimson", linewidth=0.5)
    # 标注主要频率
    for f_val, m_val in zip(dom_freqs[:5], dom_mags[:5]):
        ax3.annotate(
            f"{f_val:.1f} Hz",
            xy=(f_val, m_val),
            xytext=(f_val + 10, m_val * 1.1),
            fontsize=8,
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
            color="black",
        )
    ax3.set_title("频谱 (线性坐标)")
    ax3.set_xlabel("频率 (Hz)")
    ax3.set_ylabel("归一化幅值")
    ax3.grid(True, alpha=0.3)

    # --- 子图4: 频谱 (对数坐标) ---
    ax4 = fig.add_subplot(3, 2, 4)
    # 避免 log(0)，给 magnitude 加一个极小值
    magnitude_db = 20 * np.log10(magnitude + 1e-15)
    ax4.plot(freqs, magnitude_db, color="darkgreen", linewidth=0.5)
    ax4.set_title("频谱 (dB 标度)")
    ax4.set_xlabel("频率 (Hz)")
    ax4.set_ylabel("幅值 (dB)")
    ax4.grid(True, alpha=0.3)

    # --- 子图5: 频谱低频区域放大 (0 ~ 100 Hz) ---
    ax5 = fig.add_subplot(3, 2, 5)
    low_freq_mask = freqs <= 100
    ax5.plot(freqs[low_freq_mask], magnitude[low_freq_mask],
             color="darkorange", linewidth=1.0)
    # 标注低频区域的主要频率
    for f_val, m_val in zip(dom_freqs, dom_mags):
        if f_val <= 100:
            ax5.axvline(x=f_val, color="red", linestyle="--", alpha=0.5, linewidth=0.8)
            ax5.text(f_val + 0.5, m_val * 0.9, f"{f_val:.1f} Hz",
                     fontsize=8, color="red")
    ax5.set_title("低频区域放大 (0 ~ 100 Hz)")
    ax5.set_xlabel("频率 (Hz)")
    ax5.set_ylabel("归一化幅值")
    ax5.grid(True, alpha=0.3)

    # --- 子图6: 频谱高频区域 (100 ~ 500 Hz) ---
    ax6 = fig.add_subplot(3, 2, 6)
    high_freq_mask = (freqs >= 100) & (freqs <= 500)
    ax6.plot(freqs[high_freq_mask], magnitude[high_freq_mask],
             color="purple", linewidth=1.0)
    for f_val, m_val in zip(dom_freqs, dom_mags):
        if 100 <= f_val <= 500:
            ax6.axvline(x=f_val, color="red", linestyle="--", alpha=0.5, linewidth=0.8)
            ax6.text(f_val + 1, m_val * 0.9, f"{f_val:.1f} Hz",
                     fontsize=8, color="red")
    ax6.set_title("高频区域 (100 ~ 500 Hz)")
    ax6.set_xlabel("频率 (Hz)")
    ax6.set_ylabel("归一化幅值")
    ax6.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_figures:
        output_path = OUTPUT_DIR / "a3_signal_analysis.png"
        fig.savefig(output_path, dpi=figure_dpi, bbox_inches="tight")
        print(f"[保存] 综合分析图已保存至: {output_path}")

    plt.show()

    # ----------------------------------------------------------
    # 5. 额外可视化: 时频分析 (短时傅里叶变换频谱图)
    # ----------------------------------------------------------
    print("\n[可视化] 生成时频分析频谱图 (Spectrogram)...")

    fig2, ax_spec = plt.subplots(figsize=(14, 6))
    # 使用 matplotlib 内置的 specgram 函数
    Pxx, freqs_spec, bins_spec, im = ax_spec.specgram(
        voltage_arr,
        NFFT=1024,      # 每段 FFT 的点数
        Fs=fs,           # 采样率
        noverlap=512,    # 重叠点数
        cmap="viridis",
    )
    ax_spec.set_title("时频分析 (Spectrogram)", fontsize=14, fontweight="bold")
    ax_spec.set_xlabel("时间 (s)")
    ax_spec.set_ylabel("频率 (Hz)")
    ax_spec.set_ylim(0, 500)  # 仅显示 0~500 Hz
    fig2.colorbar(im, ax=ax_spec, label="功率谱密度 (dB/Hz)")

    plt.tight_layout()

    if save_figures:
        output_path2 = OUTPUT_DIR / "a3_spectrogram.png"
        fig2.savefig(output_path2, dpi=figure_dpi, bbox_inches="tight")
        print(f"[保存] 时频分析图已保存至: {output_path2}")

    plt.show()


if __name__ == "__main__":
    main()
