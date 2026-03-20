# -*- coding: utf-8 -*-
"""
信号分析模块

提供基于 FFT 的频域分析工具函数，用于 A.3 题的波形数据分析。

作者: kyksj-1
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple


def load_waveform(filepath: str) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    加载示波器波形数据文件。

    文件格式: 两列制表符分隔，第一列时间 (s)，第二列电压 (V)。

    Parameters
    ----------
    filepath : str
        数据文件路径。

    Returns
    -------
    Tuple[NDArray[np.float64], NDArray[np.float64]]
        (time, voltage) 元组。
        - time: 采样时间轴 (s)
        - voltage: 对应电压值 (V)
    """
    data: NDArray[np.float64] = np.loadtxt(filepath, dtype=np.float64)
    time: NDArray[np.float64] = data[:, 0]
    voltage: NDArray[np.float64] = data[:, 1]
    return time, voltage


def compute_frequency_spectrum(
    time: NDArray[np.float64],
    voltage: NDArray[np.float64],
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.complex128]]:
    """
    使用 numpy.fft 计算信号的频谱。

    仅返回正频率部分（单边谱），因为实信号的频谱关于零频对称。

    Parameters
    ----------
    time : NDArray[np.float64]
        等间距采样时间轴。
    voltage : NDArray[np.float64]
        对应的电压信号。

    Returns
    -------
    Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.complex128]]
        (freqs, magnitude, spectrum) 元组:
        - freqs: 正频率轴 (Hz)
        - magnitude: 幅值谱 |X(f)| / N（归一化后的幅值）
        - spectrum: 完整的复数频谱（正频率部分）
    """
    N: int = len(voltage)

    # 从时间轴推算采样间隔
    dt: float = time[1] - time[0]
    # 采样率
    fs: float = 1.0 / dt

    # FFT 计算（使用标准库）
    full_spectrum: NDArray[np.complex128] = np.fft.fft(voltage)
    full_freqs: NDArray[np.float64] = np.fft.fftfreq(N, d=dt)

    # 取正频率部分（单边谱）
    # rfftfreq 给出 [0, fs/2] 的频率点
    positive_mask: NDArray[np.bool_] = full_freqs >= 0
    freqs: NDArray[np.float64] = full_freqs[positive_mask]
    spectrum: NDArray[np.complex128] = full_spectrum[positive_mask]

    # 幅值谱：归一化（除以 N），乘以 2 补偿单边谱（DC 分量除外）
    magnitude: NDArray[np.float64] = np.abs(spectrum) / N
    magnitude[1:] *= 2  # 非 DC 分量乘 2

    return freqs, magnitude, spectrum


def find_dominant_frequencies(
    freqs: NDArray[np.float64],
    magnitude: NDArray[np.float64],
    threshold_ratio: float = 0.05,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    找出频谱中的主要频率分量。

    通过设定阈值（相对于最大幅值的比例）来筛选显著的频率峰。

    Parameters
    ----------
    freqs : NDArray[np.float64]
        频率轴 (Hz)。
    magnitude : NDArray[np.float64]
        幅值谱。
    threshold_ratio : float
        阈值比例，低于 max(magnitude) * threshold_ratio 的峰将被忽略。
        默认 0.05 (5%)。

    Returns
    -------
    Tuple[NDArray[np.float64], NDArray[np.float64]]
        (dominant_freqs, dominant_magnitudes):
        - dominant_freqs: 主要频率 (Hz)，按幅值降序排列
        - dominant_magnitudes: 对应的幅值
    """
    threshold: float = np.max(magnitude) * threshold_ratio

    # 简单峰检测: 幅值大于阈值且大于左右邻居
    peak_indices = []
    for i in range(1, len(magnitude) - 1):
        if (
            magnitude[i] > threshold
            and magnitude[i] > magnitude[i - 1]
            and magnitude[i] > magnitude[i + 1]
        ):
            peak_indices.append(i)

    if not peak_indices:
        return np.array([]), np.array([])

    peak_indices_arr: NDArray[np.intp] = np.array(peak_indices)

    # 按幅值降序排列
    sorted_order = np.argsort(magnitude[peak_indices_arr])[::-1]
    dominant_freqs: NDArray[np.float64] = freqs[peak_indices_arr[sorted_order]]
    dominant_magnitudes: NDArray[np.float64] = magnitude[peak_indices_arr[sorted_order]]

    return dominant_freqs, dominant_magnitudes
