# -*- coding: utf-8 -*-
"""
快速傅里叶变换 (FFT) 模块 — Base-2 Cooley-Tukey 算法

实现基于分治递归的 Radix-2 FFT，将 O(N^2) 的 DFT 降低到 O(N log N)。

核心思想 (Cooley-Tukey, 1965):
    将长度为 N 的 DFT 拆分为两个长度为 N/2 的 DFT：
    - 偶数下标子序列: x[0], x[2], x[4], ...
    - 奇数下标子序列: x[1], x[3], x[5], ...

    X[k] = E[k] + W_N^k · O[k]           (k = 0, ..., N/2 - 1)
    X[k + N/2] = E[k] - W_N^k · O[k]     (蝶形运算)

    其中 W_N^k = exp(-2πi·k / N) 为旋转因子 (twiddle factor)。

作者: kyksj-1
"""

import numpy as np
from numpy.typing import NDArray


def _is_power_of_two(n: int) -> bool:
    """
    检查整数 n 是否为 2 的幂次。

    利用位运算: 2 的幂次在二进制中只有一个 1，
    因此 n & (n - 1) == 0。

    Parameters
    ----------
    n : int
        待检查的正整数。

    Returns
    -------
    bool
        True 如果 n 是 2 的幂次，否则 False。
    """
    return n > 0 and (n & (n - 1)) == 0


def fft(x: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """
    Base-2 快速傅里叶变换 (Cooley-Tukey 递归实现)。

    复杂度: O(N log N)
    要求: 输入数组长度必须是 2 的幂次。

    Parameters
    ----------
    x : NDArray[np.complex128]
        输入的一维复数数组，长度必须为 2 的幂次。

    Returns
    -------
    NDArray[np.complex128]
        FFT 变换结果，长度为 N 的复数数组。

    Raises
    ------
    TypeError
        输入不是 numpy 数组时抛出。
    ValueError
        输入不是一维数组或长度不是 2 的幂次时抛出。
    """
    # --- 输入校验 ---
    if not isinstance(x, np.ndarray):
        raise TypeError(f"输入必须是 numpy 数组，实际类型: {type(x)}")
    if x.ndim != 1:
        raise ValueError(f"输入必须是一维数组，实际维度: {x.ndim}")

    x = x.astype(np.complex128, copy=False)
    N: int = len(x)

    if N == 0:
        return np.array([], dtype=np.complex128)

    if not _is_power_of_two(N):
        raise ValueError(
            f"输入数组长度必须是 2 的幂次，实际长度: {N}。"
            f"提示: 可以对输入进行零填充 (zero-padding) 至最近的 2 的幂次。"
        )

    return _fft_recursive(x)


def _fft_recursive(x: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """
    FFT 的递归核心实现。

    递归终止条件: 当数组长度为 1 时，DFT 就是它本身。
    递归步骤: 将数组按奇偶下标拆分，分别做 FFT，然后通过蝶形运算合并。

    Parameters
    ----------
    x : NDArray[np.complex128]
        输入数组，长度为 2 的幂次（由调用方保证）。

    Returns
    -------
    NDArray[np.complex128]
        FFT 变换结果。
    """
    N: int = len(x)

    # --- 递归终止条件 ---
    if N == 1:
        return x.copy()

    # --- 分治: 按奇偶下标拆分 ---
    x_even: NDArray[np.complex128] = x[0::2]  # 偶数下标: x[0], x[2], x[4], ...
    x_odd: NDArray[np.complex128] = x[1::2]   # 奇数下标: x[1], x[3], x[5], ...

    # 递归计算子问题的 FFT
    E: NDArray[np.complex128] = _fft_recursive(x_even)  # 偶数部分的 DFT
    O: NDArray[np.complex128] = _fft_recursive(x_odd)   # 奇数部分的 DFT

    # --- 蝶形运算 (Butterfly Operation) ---
    # 旋转因子: W_N^k = exp(-2πi·k / N), k = 0, 1, ..., N/2-1
    half_N: int = N // 2
    k: NDArray[np.float64] = np.arange(half_N)
    twiddle: NDArray[np.complex128] = np.exp(-2j * np.pi * k / N)

    # 合并: 利用 DFT 的周期性和对称性
    # X[k]       = E[k] + W_N^k · O[k]
    # X[k + N/2] = E[k] - W_N^k · O[k]
    twiddle_times_odd: NDArray[np.complex128] = twiddle * O

    X: NDArray[np.complex128] = np.empty(N, dtype=np.complex128)
    X[:half_N] = E + twiddle_times_odd
    X[half_N:] = E - twiddle_times_odd

    return X


def ifft(X: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """
    Base-2 逆快速傅里叶变换 (IFFT)。

    利用性质: IFFT(X) = conj(FFT(conj(X))) / N

    Parameters
    ----------
    X : NDArray[np.complex128]
        频域信号，长度为 2 的幂次。

    Returns
    -------
    NDArray[np.complex128]
        时域信号。
    """
    N: int = len(X)
    # IFFT 可以通过 FFT 实现: x = conj(FFT(conj(X))) / N
    return np.conj(fft(np.conj(X))) / N
