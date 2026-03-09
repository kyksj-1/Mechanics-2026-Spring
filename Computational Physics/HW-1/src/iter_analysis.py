'''
迭代映射数值分析：分岔图、Lyapunov指数、逃逸图
所有计算对 mu 向量化。
'''
import numpy as np


def compute_bifurcation(iter_func, mu_array, x0=0.1,
                         n_total=1000, n_last=200, bound=1e6):
    """计算分岔图数据，返回 (mu_points, x_points)"""
    x = np.full_like(mu_array, x0, dtype=float)

    # 丢弃瞬态
    for _ in range(n_total - n_last):
        x = iter_func(x, mu_array)
        x[np.abs(x) > bound] = np.nan

    # 采样吸引子
    mu_list, x_list = [], []
    for _ in range(n_last):
        x = iter_func(x, mu_array)
        valid = np.isfinite(x) & (np.abs(x) <= bound)
        mu_list.append(mu_array[valid])
        x_list.append(x[valid])
        x[~valid] = np.nan

    return np.concatenate(mu_list), np.concatenate(x_list)


def compute_lyapunov(iter_func, deriv_func, mu_array, x0=0.1,
                      n_total=1000, n_transient=300, bound=1e6):
    """计算 Lyapunov 指数 lambda(mu) = (1/N) sum ln|f'(x_i)|"""
    x = np.full_like(mu_array, x0, dtype=float)

    for _ in range(n_transient):
        x = iter_func(x, mu_array)
        x[np.abs(x) > bound] = np.nan

    n_sum = n_total - n_transient
    lyap = np.zeros_like(mu_array, dtype=float)
    lyap[~np.isfinite(x)] = np.nan

    for _ in range(n_sum):
        df = deriv_func(x, mu_array)
        with np.errstate(divide='ignore', invalid='ignore'):
            lyap += np.where(np.isfinite(x), np.log(np.abs(df) + 1e-30), 0)
        x = iter_func(x, mu_array)
        escaped = np.isfinite(lyap) & (~np.isfinite(x) | (np.abs(x) > bound))
        lyap[escaped] = np.nan
        x[escaped] = np.nan

    return lyap / n_sum


def compute_escape_map(iter_func, mu_array, x0_array,
                        n_total=500, bound=1e6):
    """计算 (mu, x0) 空间逃逸图，返回 bounded[i_mu, j_x0]"""
    n_mu = len(mu_array)
    bounded = np.ones((n_mu, len(x0_array)), dtype=bool)

    for j, x0 in enumerate(x0_array):
        x = np.full(n_mu, x0, dtype=float)
        for _ in range(n_total):
            x = iter_func(x, mu_array)
            escaped = np.abs(x) > bound
            bounded[escaped, j] = False
            x[escaped] = np.nan

    return bounded
