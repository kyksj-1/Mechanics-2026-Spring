# -*- coding: utf-8 -*-
"""
牛顿迭代法模块 — 复平面上的非线性方程求解与 Newton 分形

实现通用的牛顿迭代法框架，并针对 f(z) = z^3 - 1 提供高性能的
向量化分形计算。

数学背景:
    牛顿迭代公式: z_{n+1} = z_n - f(z_n) / f'(z_n)
    对于 f(z) = z^3 - 1:
        f'(z) = 3z^2
        z_{n+1} = z_n - (z_n^3 - 1) / (3 z_n^2)
                = (2 z_n^3 + 1) / (3 z_n^2)

    三个根 (三次单位根):
        r_0 = 1
        r_1 = exp(2πi/3) = -1/2 + i√3/2
        r_2 = exp(4πi/3) = -1/2 - i√3/2

作者: kyksj-1
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, Callable, Optional
import logging

logger = logging.getLogger(__name__)


# ============================================================
# z^3 - 1 的三个根 (三次单位根)
# ============================================================
ROOTS_Z3: NDArray[np.complex128] = np.array([
    1.0 + 0.0j,                                    # r_0 = 1
    -0.5 + np.sqrt(3) / 2 * 1j,                    # r_1 = exp(2πi/3)
    -0.5 - np.sqrt(3) / 2 * 1j,                    # r_2 = exp(4πi/3)
], dtype=np.complex128)


def newton_step_z3(z: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """
    f(z) = z^3 - 1 的单步牛顿迭代（向量化）。

    推导:
        z_{n+1} = z_n - f(z_n)/f'(z_n)
                = z_n - (z_n^3 - 1) / (3 z_n^2)

    为避免除零，对 |z| 极小的点做保护处理。

    Parameters
    ----------
    z : NDArray[np.complex128]
        当前迭代点（可以是任意形状的数组）。

    Returns
    -------
    NDArray[np.complex128]
        迭代一步后的 z 值。
    """
    z_sq: NDArray[np.complex128] = z * z           # z^2
    denom: NDArray[np.complex128] = 3.0 * z_sq     # f'(z) = 3z^2

    # 对 |f'(z)| 极小的点进行保护，防止除零导致 NaN/Inf
    safe_mask: NDArray[np.bool_] = np.abs(denom) > 1e-15
    # 对不安全的点，保持原值不动（即跳过迭代）
    result: NDArray[np.complex128] = z.copy()
    result[safe_mask] = z[safe_mask] - (z[safe_mask] ** 3 - 1.0) / denom[safe_mask]

    return result


def compute_newton_fractal(
    center: Tuple[float, float],
    half_width: float,
    resolution: float,
    max_iter: int = 100,
    tol: float = 1e-8,
    roots: Optional[NDArray[np.complex128]] = None,
) -> Tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.float64], NDArray[np.float64]]:
    """
    计算 Newton 分形: 对网格上每个点执行牛顿迭代，记录收敛到哪个根。

    Parameters
    ----------
    center : Tuple[float, float]
        区域中心坐标 (cx, cy)。
    half_width : float
        正方形区域的半宽。
    resolution : float
        网格分辨率（像素间距）。
    max_iter : int
        最大迭代次数，默认 100。
    tol : float
        收敛判据: |z_n - root| < tol 时认为收敛。
    roots : Optional[NDArray[np.complex128]]
        方程的根。默认为 z^3-1 的三个根。

    Returns
    -------
    Tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.float64], NDArray[np.float64]]
        (root_index, iter_count, x_axis, y_axis):
        - root_index: 形状 (ny, nx) 的整数数组，值为收敛到的根的编号
          (0, 1, 2)，未收敛的点标记为 -1
        - iter_count: 形状 (ny, nx) 的整数数组，收敛所需的迭代次数
        - x_axis: x 方向的坐标数组
        - y_axis: y 方向的坐标数组
    """
    if roots is None:
        roots = ROOTS_Z3

    cx, cy = center

    # --- 构建网格 ---
    x_axis: NDArray[np.float64] = np.arange(
        cx - half_width, cx + half_width + resolution * 0.5, resolution
    )
    y_axis: NDArray[np.float64] = np.arange(
        cy - half_width, cy + half_width + resolution * 0.5, resolution
    )

    nx, ny = len(x_axis), len(y_axis)
    logger.info(f"网格大小: {nx} x {ny} = {nx * ny} 个点")

    # 构造复数网格 z = x + iy
    # meshgrid 返回 (ny, nx) 形状（行对应 y，列对应 x）
    X, Y = np.meshgrid(x_axis, y_axis)
    Z: NDArray[np.complex128] = (X + 1j * Y).astype(np.complex128)

    # --- 初始化结果数组 ---
    # root_index[i,j] = 收敛到的根编号，-1 表示未收敛
    root_index: NDArray[np.int32] = np.full(Z.shape, -1, dtype=np.int32)
    # iter_count[i,j] = 收敛所需迭代次数
    iter_count: NDArray[np.int32] = np.full(Z.shape, max_iter, dtype=np.int32)

    # 活跃点掩码: True 表示该点尚未收敛
    active: NDArray[np.bool_] = np.ones(Z.shape, dtype=bool)

    # --- 牛顿迭代 ---
    for iteration in range(max_iter):
        # 仅对活跃点执行迭代
        Z[active] = newton_step_z3(Z[active])

        # 检查收敛: 对每个根计算距离
        for root_idx, root in enumerate(roots):
            # 距离小于 tol 的点视为收敛到该根
            converged: NDArray[np.bool_] = active & (np.abs(Z - root) < tol)
            root_index[converged] = root_idx
            iter_count[converged] = iteration + 1
            active[converged] = False

        # 若所有点均已收敛，提前退出
        if not np.any(active):
            logger.info(f"所有点在第 {iteration + 1} 次迭代后收敛")
            break

    # 统计未收敛点数
    n_unconverged: int = int(np.sum(active))
    if n_unconverged > 0:
        logger.warning(f"{n_unconverged} 个点在 {max_iter} 次迭代后未收敛")

    return root_index, iter_count, x_axis, y_axis


def newton_iterate_general(
    f: Callable[[complex], complex],
    f_prime: Callable[[complex], complex],
    z0: complex,
    max_iter: int = 100,
    tol: float = 1e-10,
) -> Tuple[complex, int, bool]:
    """
    通用牛顿迭代法（标量版本，用于单点演示/验证）。

    z_{n+1} = z_n - f(z_n) / f'(z_n)

    Parameters
    ----------
    f : Callable[[complex], complex]
        目标函数。
    f_prime : Callable[[complex], complex]
        目标函数的导数。
    z0 : complex
        初始猜测值。
    max_iter : int
        最大迭代次数。
    tol : float
        收敛容差。

    Returns
    -------
    Tuple[complex, int, bool]
        (z_final, n_iterations, converged):
        - z_final: 最终迭代结果
        - n_iterations: 实际迭代次数
        - converged: 是否收敛
    """
    z: complex = z0
    for i in range(max_iter):
        fz: complex = f(z)
        fpz: complex = f_prime(z)

        if abs(fpz) < 1e-15:
            # 导数过小，无法继续迭代
            return z, i, False

        z_new: complex = z - fz / fpz

        if abs(z_new - z) < tol:
            return z_new, i + 1, True

        z = z_new

    return z, max_iter, False
