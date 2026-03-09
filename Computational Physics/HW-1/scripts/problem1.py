'''
题1: 虫口模型与 Logistic 映射
Q1: x_{n+1} = 1 - mu * x_n^2,       mu in (0,2), x0 in (-1,1)
Q2: x_{n+1} = cos(x_n) - mu * x_n^2, 合理选择区间
Q3: 两种映射对比
'''
import numpy as np
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from src.iterative_equation import (
    linear_iter, linear_deriv,
    cosine_iter, cosine_deriv,
)
from src.iter_analysis import compute_bifurcation, compute_lyapunov
from src.iter_plot import plot_combined, plot_timeseries, plot_bifurcation

# ---- 输出目录 ----
OUTDIR = ROOT / 'plots' / 'asset' / 'q1'
OUTDIR.mkdir(parents=True, exist_ok=True)

# ---- 计算参数 ----
N_TOTAL = 2000       # 总迭代步数
N_LAST = 300         # 分岔图采样最后 300 步
N_TRANSIENT = 500    # Lyapunov 丢弃瞬态步数
BOUND = 1e6          # 发散阈值


def analyze_map(tag, label, iter_func, deriv_func, mu_range,
                x0_bif=0.1, ts_mus=None, ts_x0s=None):
    """
    通用分析管线：分岔图 + Lyapunov、逃逸图、时间序列
    tag: 文件名前缀 (e.g. 'linear')
    label: 图标题中的方程表达式
    返回: (mu_bif, x_bif) 供后续对比使用
    """
    print(f'\n{"="*50}')
    print(f'  {tag}')
    print(f'{"="*50}')

    mu_fine = np.linspace(mu_range[0], mu_range[1], 2000)

    # 1. 分岔图
    print('  [1/4] 分岔图...')
    mu_bif, x_bif = compute_bifurcation(
        iter_func, mu_fine, x0=x0_bif,
        n_total=N_TOTAL, n_last=N_LAST, bound=BOUND
    )

    # 2. Lyapunov 指数
    print('  [2/4] Lyapunov 指数...')
    lyap = compute_lyapunov(
        iter_func, deriv_func, mu_fine, x0=x0_bif,
        n_total=N_TOTAL, n_transient=N_TRANSIENT, bound=BOUND
    )

    # 组合图：分岔 + Lyapunov
    fig = plot_combined(mu_bif, x_bif, mu_fine, lyap, suptitle=label)
    path = OUTDIR / f'{tag}_bifurcation_lyapunov.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'    -> {path}')
    plt.close(fig)

    # 3. 初值依赖性：多 x0 分岔图叠加
    print('  [3/4] 初值依赖性...')
    x0_test = [0.1, 0.5, -0.3]
    colors = ['C0', 'C1', 'C2']
    fig, ax = plt.subplots(figsize=(12, 7))
    for x0_val, c in zip(x0_test, colors):
        mu_b, x_b = compute_bifurcation(
            iter_func, mu_fine, x0=x0_val,
            n_total=N_TOTAL, n_last=N_LAST, bound=BOUND
        )
        ax.scatter(mu_b, x_b, s=0.02, c=c, alpha=0.4,
                   edgecolors='none', rasterized=True,
                   label=f'$x_0={x0_val}$')
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'$x^*$')
    ax.set_title(f'Initial condition dependence ({tag})')
    ax.legend(markerscale=40, loc='upper left')
    path = OUTDIR / f'{tag}_x0_dependence.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'    -> {path}')
    plt.close(fig)

    # 4. 时间序列
    if ts_mus is not None and ts_x0s is not None:
        print('  [4/4] 时间序列...')
        fig = plot_timeseries(
            iter_func, ts_mus, ts_x0s,
            n_iter=200, bound=BOUND, suptitle=label
        )
        path = OUTDIR / f'{tag}_timeseries.png'
        fig.savefig(path, dpi=200, bbox_inches='tight')
        print(f'    -> {path}')
        plt.close(fig)

    return mu_bif, x_bif


if __name__ == '__main__':
    # ===== Q1: Linear 映射 =====
    bif_l = analyze_map(
        tag='linear',
        label=r'$x_{n+1}=1-\mu x_n^2$',
        iter_func=linear_iter,
        deriv_func=linear_deriv,
        mu_range=(0.01, 2.0),
        ts_mus=[0.5, 0.85, 1.25, 1.8],
        ts_x0s=[-0.5, 0.1, 0.8],
    )

    # ===== Q2: Cosine 映射 =====
    # cosine 首次分岔 ~ mu=0.28，远早于 linear 的 0.75
    # mu 取 (0, 1.2) 覆盖完整的倍周期→混沌→逃逸过程
    bif_c = analyze_map(
        tag='cosine',
        label=r'$x_{n+1}=\cos x_n - \mu x_n^2$',
        iter_func=cosine_iter,
        deriv_func=cosine_deriv,
        mu_range=(0.01, 1.2),
        ts_mus=[0.1, 0.35, 0.55, 0.8],
        ts_x0s=[-0.5, 0.1, 0.8],
    )

    # ===== Q3: 并排对比分岔图 =====
    print('\n[Q3] 对比...')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    plot_bifurcation(bif_l[0], bif_l[1],
                     title=r'$x_{n+1}=1-\mu x_n^2$', ax=ax1)
    plot_bifurcation(bif_c[0], bif_c[1],
                     title=r'$x_{n+1}=\cos x_n - \mu x_n^2$', ax=ax2)
    fig.suptitle('Bifurcation Diagram Comparison', fontsize=14)
    fig.tight_layout()
    path = OUTDIR / 'comparison_bifurcation.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  -> {path}')
    plt.close(fig)

    print('\n全部完成!')
