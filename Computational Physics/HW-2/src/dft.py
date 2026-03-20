# -*- coding: utf-8 -*-
"""
离散傅里叶变换 (DFT) 模块

实现标准的 O(N^2) 离散傅里叶变换，严格按照定义：
    X[k] = Σ_{n=0}^{N-1} x[n] · exp(-2πi·n·k / N),  k = 0, 1, ..., N-1

作者: kyksj-1
"""

import numpy as np
from numpy.typing import NDArray


def dft(x: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """
    计算一维离散傅里叶变换 (DFT)。

    严格按照 DFT 定义式，使用向量化矩阵乘法实现。
    复杂度: O(N^2)

    Parameters
    ----------
    x : NDArray[np.complex128]
        输入的一维复数数组，长度为 N。

    Returns
    -------
    NDArray[np.complex128]
        DFT 变换结果，长度为 N 的复数数组。

    Raises
    ------
    TypeError
        输入不是 numpy 数组时抛出。
    ValueError
        输入不是一维数组时抛出。
    """
    # --- 输入校验 ---
    if not isinstance(x, np.ndarray):
        raise TypeError(f"输入必须是 numpy 数组，实际类型: {type(x)}")
    if x.ndim != 1:
        raise ValueError(f"输入必须是一维数组，实际维度: {x.ndim}")

    # 确保输入为 complex128 类型
    x = x.astype(np.complex128, copy=False)
    N: int = len(x)

    if N == 0:
        return np.array([], dtype=np.complex128)

    # --- 构造 DFT 矩阵 ---
    # n, k 分别为行索引和列索引（均从 0 到 N-1）
    n: NDArray[np.float64] = np.arange(N)
    k: NDArray[np.float64] = np.arange(N)

    # DFT 矩阵元素: W[k, n] = exp(-2πi·n·k / N)
    # 利用外积构造 N×N 矩阵，其中 (k, n) 元素为 k*n
    exponent_matrix: NDArray[np.complex128] = np.exp(
        -2j * np.pi * np.outer(k, n) / N
    )

    # --- 矩阵-向量乘法得到 DFT 结果 ---
    # X = W @ x，即 X[k] = Σ_n W[k,n] · x[n]
    X: NDArray[np.complex128] = exponent_matrix @ x

    return X


def idft(X: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """
    计算一维逆离散傅里叶变换 (IDFT)。

    定义: x[n] = (1/N) Σ_{k=0}^{N-1} X[k] · exp(2πi·n·k / N)

    Parameters
    ----------
    X : NDArray[np.complex128]
        频域信号，长度为 N 的复数数组。

    Returns
    -------
    NDArray[np.complex128]
        时域信号，长度为 N 的复数数组。
    """
    if not isinstance(X, np.ndarray):
        raise TypeError(f"输入必须是 numpy 数组，实际类型: {type(X)}")
    if X.ndim != 1:
        raise ValueError(f"输入必须是一维数组，实际维度: {X.ndim}")

    X = X.astype(np.complex128, copy=False)
    N: int = len(X)

    if N == 0:
        return np.array([], dtype=np.complex128)

    n: NDArray[np.float64] = np.arange(N)
    k: NDArray[np.float64] = np.arange(N)

    # IDFT 矩阵: W_inv[n, k] = exp(+2πi·n·k / N)
    # 注意指数符号为正
    exponent_matrix: NDArray[np.complex128] = np.exp(
        2j * np.pi * np.outer(n, k) / N
    )

    x: NDArray[np.complex128] = (exponent_matrix @ X) / N
    return x
