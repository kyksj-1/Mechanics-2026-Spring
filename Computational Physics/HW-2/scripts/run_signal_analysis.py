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
    ax1.set_title("完整时域波形 (40s)")
    ax1.set_xlabel("时间 (s)")
    ax1.set_ylabel("电压 (V)")
    ax1.grid(True, alpha=0.3)

    # --- 子图2: 时域波形局部放大 (前 5s, 约 5 个完整周期) ---
    ax2 = fig.add_subplot(3, 2, 2)
    mask_zoom = time_arr <= 5.0
    ax2.plot(time_arr[mask_zoom], voltage_arr[mask_zoom],
             color="steelblue", linewidth=0.5)
    ax2.set_title("时域波形 (前 5s, 约 5 个基频周期)")
    ax2.set_xlabel("时间 (s)")
    ax2.set_ylabel("电压 (V)")
    ax2.grid(True, alpha=0.3)

    # --- 子图3: 频谱 — 信号区域 (0 ~ 15 Hz) ---
    ax3 = fig.add_subplot(3, 2, 3)
    sig_mask = freqs <= 15
    ax3.plot(freqs[sig_mask], magnitude[sig_mask], color="crimson", linewidth=1.2)
    # 用茎叶图 (stem) 标注主要频率峰，使其更醒目
    for f_val, m_val in zip(dom_freqs, dom_mags):
        if f_val <= 15:
            ax3.plot(f_val, m_val, "o", color="darkred", markersize=7, zorder=5)
            ax3.annotate(
                f"{f_val:.0f} Hz\n(A={m_val:.3f})",
                xy=(f_val, m_val),
                xytext=(f_val + 0.8, m_val - 0.06),
                fontsize=9, fontweight="bold", color="darkred",
                arrowprops=dict(arrowstyle="->", color="darkred", lw=1.0),
            )
    ax3.set_title("频谱 — 信号区域 (0~15 Hz)")
    ax3.set_xlabel("频率 (Hz)")
    ax3.set_ylabel("归一化幅值")
    ax3.set_ylim(-0.02, 1.1)
    ax3.grid(True, alpha=0.3)

    # --- 子图4: 频谱 (dB 标度, 0~500 Hz) ---
    ax4 = fig.add_subplot(3, 2, 4)
    magnitude_db = 20 * np.log10(magnitude + 1e-15)
    ax4.plot(freqs, magnitude_db, color="darkgreen", linewidth=0.5)
    # 标注主要峰
    for f_val, m_val in zip(dom_freqs, dom_mags):
        db_val = 20 * np.log10(m_val + 1e-15)
        ax4.plot(f_val, db_val, "v", color="red", markersize=6, zorder=5)
    ax4.set_title("频谱 (dB 标度, 全频段)")
    ax4.set_xlabel("频率 (Hz)")
    ax4.set_ylabel("幅值 (dB)")
    ax4.grid(True, alpha=0.3)

    # --- 子图5: 与方波傅里叶级数的理论比较 ---
    ax5 = fig.add_subplot(3, 2, 5)
    # 理论方波谐波: A_n = (4/pi) * (1/n), n = 1, 3, 5, 7, ...
    theory_n = np.array([1, 3, 5, 7, 9, 11, 13])
    theory_amp = 1.0 / theory_n  # 归一化到基频=1
    # 实测幅值 (归一化到基频=1)
    measured_amp = dom_mags / dom_mags[0] if len(dom_mags) > 0 else np.array([])
    measured_n = dom_freqs / dom_freqs[0] if len(dom_freqs) > 0 else np.array([])

    bar_width = 0.35
    x_pos = np.arange(len(theory_n))
    ax5.bar(x_pos - bar_width/2, theory_amp[:len(theory_n)], bar_width,
            color="steelblue", alpha=0.8, label="理论 (1/n)")
    # 只画有实测数据的谐波
    n_measured = min(len(measured_amp), len(theory_n))
    ax5.bar(x_pos[:n_measured] + bar_width/2, measured_amp[:n_measured], bar_width,
            color="coral", alpha=0.8, label="实测")
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels([f"{int(n)}f0" for n in theory_n])
    ax5.set_title("谐波幅值: 实测 vs 方波理论 (1/n)")
    ax5.set_xlabel("谐波次数")
    ax5.set_ylabel("归一化幅值 (基频=1)")
    ax5.legend(fontsize=10)
    ax5.grid(True, alpha=0.3, axis="y")

    # --- 子图6: 噪声本底分析 (10 ~ 500 Hz) ---
    ax6 = fig.add_subplot(3, 2, 6)
    noise_mask = freqs >= 10
    ax6.plot(freqs[noise_mask], magnitude[noise_mask],
             color="purple", linewidth=0.5, alpha=0.7)
    noise_floor = np.median(magnitude[noise_mask])
    ax6.axhline(y=noise_floor, color="red", linestyle="--", linewidth=1.5,
                label=f"噪声中位值 = {noise_floor:.4f}")
    ax6.set_title("噪声本底 (10~500 Hz, 信号峰之外)")
    ax6.set_xlabel("频率 (Hz)")
    ax6.set_ylabel("归一化幅值")
    ax6.legend(fontsize=10)
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
        NFFT=1024,      # 每段 FFT 的点数（频率分辨率 ≈ fs/NFFT ≈ 0.98 Hz）
        Fs=fs,           # 采样率
        noverlap=512,    # 重叠点数（50% 重叠）
        cmap="viridis",
    )
    ax_spec.set_title("时频分析 (Spectrogram)", fontsize=14, fontweight="bold")
    ax_spec.set_xlabel("时间 (s)")
    ax_spec.set_ylabel("频率 (Hz)")
    # 信号分量集中在 1/3/5/7 Hz，将 y 轴聚焦到信号区域以清晰显示各谐波条带
    ax_spec.set_ylim(0, 15)
    # 在各谐波频率处添加参考虚线，便于对照
    for harmonic_freq in [1, 3, 5, 7]:
        ax_spec.axhline(
            y=harmonic_freq, color="white", linestyle="--",
            linewidth=0.8, alpha=0.7,
        )
        ax_spec.text(
            0.3, harmonic_freq + 0.3, f"{harmonic_freq} Hz",
            color="white", fontsize=9, fontweight="bold", alpha=0.9,
        )
    fig2.colorbar(im, ax=ax_spec, label="功率谱密度 (dB/Hz)")

    plt.tight_layout()

    if save_figures:
        output_path2 = OUTPUT_DIR / "a3_spectrogram.png"
        fig2.savefig(output_path2, dpi=figure_dpi, bbox_inches="tight")
        print(f"[保存] 时频分析图已保存至: {output_path2}")

    plt.show()


if __name__ == "__main__":
    main()
