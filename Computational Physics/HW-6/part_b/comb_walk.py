"""
梳子晶格随机行走矢量模拟器

接口说明
--------
simulate_comb_walk(walkers, max_steps) -> (x, y)
    - 矢量化并行模拟 N 个行走者在梳子上的随机行走
    - 返回每个行走者每步的位置

simulate_1d_walk(walkers, max_steps) -> x
    - 一维随机行走（用于对比）

设计要点
--------
- 每步对全体行走者矢量化操作
- 用布尔索引区分"在主干上"和"在分支上"的行走者
- 避免逐行走者的 Python 循环
"""

import numpy as np
from typing import Tuple


def simulate_comb_walk(
    num_walkers: int,
    max_steps: int,
    seed: int = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    矢量化梳子晶格随机行走

    规则:
      - y=0 (主干上): 等概率选 ±x, ±y 四个方向
      - y≠0 (分支上): 等概率选 ±y 两个方向

    参数:
        num_walkers : 行走者数量
        max_steps   : 每行走者的最大步数
        seed        : 随机种子

    返回:
        x : shape (num_walkers, max_steps+1) -- x 坐标轨迹
        y : shape (num_walkers, max_steps+1) -- y 坐标轨迹
    """
    if seed is not None:
        np.random.seed(seed)

    x = np.zeros((num_walkers, max_steps + 1), dtype=np.int32)
    y = np.zeros((num_walkers, max_steps + 1), dtype=np.int32)

    for step in range(max_steps):
        x_curr = x[:, step]
        y_curr = y[:, step]

        on_spine = (y_curr == 0)

        # ---- 主干上: 随机选四个方向 (0:+x, 1:-x, 2:+y, 3:-y) ----
        # ---- 分支上: 随机选两个方向 (0:+y, 1:-y) ----
        dx = np.zeros(num_walkers, dtype=np.int32)
        dy = np.zeros(num_walkers, dtype=np.int32)

        # 主干上
        if np.any(on_spine):
            spine_choice = np.random.randint(0, 4, size=on_spine.sum())
            dx[on_spine] = np.where(spine_choice == 0, 1,
                           np.where(spine_choice == 1, -1, 0))
            dy[on_spine] = np.where(spine_choice == 2, 1,
                           np.where(spine_choice == 3, -1, 0))

        # 分支上
        off_spine = ~on_spine
        if np.any(off_spine):
            branch_choice = np.random.randint(0, 2, size=off_spine.sum())
            dy[off_spine] = np.where(branch_choice == 0, 1, -1)

        x[:, step + 1] = x_curr + dx
        y[:, step + 1] = y_curr + dy

    return x, y


def simulate_1d_walk(
    num_walkers: int,
    max_steps: int,
    seed: int = None
) -> np.ndarray:
    """
    矢量化一维随机行走 (对比用)

    返回: x : shape (num_walkers, max_steps+1)
    """
    if seed is not None:
        np.random.seed(seed)

    steps = np.random.choice([-1, 1], size=(num_walkers, max_steps))
    x = np.zeros((num_walkers, max_steps + 1), dtype=np.int32)
    x[:, 1:] = np.cumsum(steps, axis=1)
    return x
