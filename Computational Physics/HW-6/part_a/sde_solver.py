"""
SDE求解器模块 —— Langevin方程的矢量化数值求解

接口说明
--------
simulate_langevin : 核心求解函数，并行模拟 N 个粒子的 Langevin 动力学
    - 输入：初始位置 q0, 初始速度 v0, 势能力函数 V_prime, 参数 params,
           时间步长 dt, 总步数 num_steps, 粒子数 num_particles
    - 输出：q, v (shape均为 (num_particles, num_steps+1))

compute_ensemble_energies : 从轨迹计算系综平均动能和势能
    - 输入：q, v, V, m
    - 输出：Ek_avg(t), V_avg(t) (1D数组)

设计要点
--------
- 所有操作均矢量化，利用 numpy 数组运算避免 Python 循环
- 使用 symplectic Euler-Maruyama 格式，保持确定性部分的辛结构
- 随机数一次性预生成，形状为 (num_particles, num_steps)
"""

import numpy as np
from typing import Callable, Tuple


def simulate_langevin(
    q0: np.ndarray,          # 初始位置, shape: (num_particles,) 
    v0: np.ndarray,          # 初始速度, shape: (num_particles,)
    V_prime: Callable,       # 势能导数 dV/dq, 支持矢量化调用
    m: float,                # 质量
    lam: float,              # 阻尼系数 λ
    k_B: float,              # 玻尔兹曼常数
    T: float,                # 温度
    dt: float,               # 时间步长
    num_steps: int,          # 总步数
    num_particles: int       # 粒子总数
) -> Tuple[np.ndarray, np.ndarray]:
    """
    矢量化模拟 Langevin 动力学 (symplectic Euler-Maruyama 格式)

    离散格式:
        v_{n+1} = v_n - γ v_n Δt - (1/m) V'(q_n) Δt + (D/m) sqrt(Δt) ζ_n
        q_{n+1} = q_n + v_{n+1} Δt

    返回: q, v -- shape (num_particles, num_steps+1)
    """
    gamma = lam / m
    D = np.sqrt(2.0 * lam * k_B * T)
    noise_std = (D / m) * np.sqrt(dt)

    # 预生成所有高斯噪声: shape (num_particles, num_steps)
    zeta = np.random.randn(num_particles, num_steps) * noise_std

    q = np.zeros((num_particles, num_steps + 1))
    v = np.zeros((num_particles, num_steps + 1))
    q[:, 0] = q0
    v[:, 0] = v0

    # 矢量化时间步进 (沿时间轴逐列计算，粒子轴已矢量化)
    q_curr = q0.copy()
    v_curr = v0.copy()
    for n in range(num_steps):
        force = -V_prime(q_curr) / m
        v_next = v_curr + (force - gamma * v_curr) * dt + zeta[:, n]
        q_curr = q_curr + v_next * dt
        v_curr = v_next
        q[:, n + 1] = q_curr
        v[:, n + 1] = v_curr

    return q, v


def compute_ensemble_energies(
    q: np.ndarray,    # shape: (num_particles, num_steps+1)
    v: np.ndarray,    # shape: (num_particles, num_steps+1)
    V: Callable,      # 势能函数 V(q)
    m: float          # 质量
) -> Tuple[np.ndarray, np.ndarray]:
    """
    计算系综平均动能和势能

    返回: Ek_avg (shape: num_steps+1), V_avg (shape: num_steps+1)
    """
    Ek = 0.5 * m * v**2
    V_q = V(q)
    return Ek.mean(axis=0), V_q.mean(axis=0)
