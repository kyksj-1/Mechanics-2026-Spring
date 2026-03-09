# src/container_parallel.py
import numpy as np
import multiprocessing as mp
from src.basic_component import GridPlaygroundContainer

def log_binning(data: np.ndarray, a: float = 1.3):
    """
    对数分箱算法：用于平滑双对数坐标下的长尾幂律分布噪声。
    参数 a 控制箱子宽度的指数增长率。
    """
    data = np.asarray(data)
    if data.size == 0:
        return np.array([]), np.array([])

    if a <= 1.0:
        raise ValueError("a must be > 1.0")

    max_val = int(np.max(data))
    if max_val < 1:
        return np.array([]), np.array([])

    edges = [1.0]
    while edges[-1] < (max_val + 1):
        edges.append(edges[-1] * a)

    edges = np.unique(np.floor(edges).astype(int))
    edges = edges[edges >= 1]
    if edges.size < 2:
        edges = np.array([1, max_val + 1], dtype=int)
    elif edges[-1] <= max_val:
        edges = np.append(edges, max_val + 1)

    hist, bin_edges = np.histogram(data, bins=edges)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    
    # 归一化为概率密度
    prob = hist / (bin_widths * len(data))
    
    # 滤除概率为0的空箱，以适配对数坐标系
    mask = prob > 0
    return bin_centers[mask], prob[mask]

def fit_powerlaw(s_centers: np.ndarray, prob: np.ndarray, fit_range: tuple = None):
    """
    在双对数坐标下对 P(s) 进行线性拟合，提取幂律指数 tau。
    纯工具函数，不依赖任何绘图或模拟逻辑。

    参数:
        s_centers: 对数分箱后的中心值
        prob: 对应的概率密度
        fit_range: (s_min, s_max) 拟合范围，默认自动选取中间段
    返回:
        tau: 幂律指数（正值，P(s) ~ s^{-tau}）
        fit_s: 拟合区间的 s 值
        fit_p: 拟合得到的 P(s) 值
    """
    log_s = np.log10(s_centers)
    log_p = np.log10(prob)

    if fit_range is not None:
        mask = (s_centers >= fit_range[0]) & (s_centers <= fit_range[1])
    else:
        # 自动选取中间段：跳过首尾的离散效应和有限尺寸截断
        n = len(log_s)
        start = max(1, int(n * 0.1))
        end = max(start + 2, int(n * 0.85))
        mask = np.zeros(n, dtype=bool)
        mask[start:end] = True

    if np.sum(mask) < 2:
        raise ValueError("拟合数据点不足，请增加样本量或调整拟合范围")

    coeffs = np.polyfit(log_s[mask], log_p[mask], 1)
    slope, intercept = coeffs
    tau = -slope

    fit_s = 10 ** log_s[mask]
    fit_p = 10 ** (slope * log_s[mask] + intercept)

    return tau, fit_s, fit_p


def worker_func(L: int, t_step_max: int, burn_in_ratio: float = 0.3):
    """进程独立工作函数，避开 Python GIL"""
    container = GridPlaygroundContainer(L)
    playground = np.zeros((L, L), dtype=int)
    
    burn_in_steps = int(t_step_max * burn_in_ratio)
    scores = []
    
    for t in range(t_step_max):
        playground, score, _ = container.one_step(playground)
        # 必须舍弃瞬态数据，只记录系统达到临界稳态后的雪崩分数
        if t > burn_in_steps and score > 0:
            scores.append(score)
            
    return L, np.array(scores)

class ContainerParallel:
    def __init__(self, Ls: list, t_step_max: int = 50000):
        self.Ls = Ls
        self.t_step_max = t_step_max
        self.results = {}
        
    def run_parallel(self):
        """利用多进程并发执行不同系统尺寸的演化"""
        print(f"开始并行计算，系统尺寸: {self.Ls}，总步数: {self.t_step_max}")
        args = [(L, self.t_step_max) for L in self.Ls]
        
        with mp.Pool(processes=len(self.Ls)) as pool:
            # starmap 展开参数列表并传入 worker_func
            parallel_results = pool.starmap(worker_func, args)
            
        for L, scores in parallel_results:
            self.results[L] = scores
        print("并行计算完成。")
