'''
迭代映射可视化：分岔图、Lyapunov指数、逃逸图、时间序列
'''
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def plot_bifurcation(mu_pts, x_pts, title='', ax=None, **kwargs):
    """散点分岔图"""
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 7))
    defaults = dict(s=0.02, c='black', alpha=0.5, edgecolors='none', rasterized=True)
    defaults.update(kwargs)
    ax.scatter(mu_pts, x_pts, **defaults)
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'$x^*$')
    if title:
        ax.set_title(title)
    return ax


def plot_lyapunov(mu_array, lyap_array, title='', ax=None):
    """Lyapunov 指数 vs mu"""
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 3))
    valid = np.isfinite(lyap_array)
    ax.plot(mu_array[valid], lyap_array[valid], lw=0.4, color='steelblue')
    ax.axhline(0, color='red', ls='--', lw=0.8, alpha=0.7)
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'$\lambda$')
    if title:
        ax.set_title(title)
    return ax


def plot_combined(mu_bif, x_bif, mu_lyap, lyap, suptitle=''):
    """分岔图(上) + Lyapunov指数(下)，共享横轴"""
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 9),
        gridspec_kw={'height_ratios': [3, 1]},
        sharex=True
    )
    plot_bifurcation(mu_bif, x_bif, ax=ax1)
    ax1.set_xlabel('')
    plot_lyapunov(mu_lyap, lyap, ax=ax2)
    if suptitle:
        fig.suptitle(suptitle, fontsize=14)
    fig.tight_layout()
    return fig


def plot_escape_map(mu_array, x0_array, bounded, title='', ax=None):
    """(mu, x0) 参数空间逃逸图，黑=有界，白=发散"""
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    extent = [mu_array[0], mu_array[-1], x0_array[0], x0_array[-1]]
    ax.imshow(bounded.T, origin='lower', aspect='auto', extent=extent,
              cmap='gray_r', interpolation='nearest', vmin=0, vmax=1)
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'$x_0$')
    if title:
        ax.set_title(title)
    return ax


def plot_timeseries(iter_func, mus, x0s, n_iter=200, bound=1e6, suptitle=''):
    """时间序列网格：行=mu，列=x0"""
    n_mu, n_x0 = len(mus), len(x0s)
    fig, axes = plt.subplots(n_mu, n_x0, figsize=(4 * n_x0, 2.5 * n_mu),
                              squeeze=False, sharex=True)
    for i, mu in enumerate(mus):
        for j, x0 in enumerate(x0s):
            ax = axes[i][j]
            traj = [x0]
            x = x0
            for _ in range(n_iter):
                x = iter_func(x, mu)
                if abs(x) > bound:
                    break
                traj.append(x)
            ax.plot(traj, lw=0.5)
            if i == 0:
                ax.set_title(f'$x_0={x0:.2f}$', fontsize=9)
            if j == 0:
                ax.set_ylabel(f'$\\mu={mu:.2f}$\n$x$', fontsize=9)
            if i == n_mu - 1:
                ax.set_xlabel('$n$', fontsize=9)
    if suptitle:
        fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout()
    return fig


def fig_save(fig=None, filepath=None,
            equ_name='linear',
            mu0s=None, x0s=None):
    '''保存图片到指定文件夹'''
    if filepath is None:
        mu_range = "mu"
        x0_range = "x0"
        if mu0s is not None:
            mu_range = f"mu[{float(np.min(mu0s)):.2f},{float(np.max(mu0s)):.2f}]"
        if x0s is not None:
            x0_range = f"x0[{float(np.min(x0s)):.2f},{float(np.max(x0s)):.2f}]"
        base_dir = Path.cwd() / "plots" / f"{equ_name}_{mu_range}_{x0_range}"
    else:
        base_dir = Path(filepath)
    base_dir.mkdir(parents=True, exist_ok=True)
    figs = [fig] if fig is not None else [plt.figure(n) for n in plt.get_fignums()]
    for f in figs:
        mu_text = None
        if hasattr(f, "_suptitle") and f._suptitle is not None:
            mu_text = f._suptitle.get_text()
        if not mu_text:
            for ax in f.get_axes():
                t = ax.get_title()
                if t.startswith("mu="):
                    mu_text = t
                    break
        mu_part = "mu" if mu_text is None else mu_text
        out_path = base_dir / f"{equ_name}_{mu_part}.png"
        f.savefig(out_path, dpi=200, bbox_inches="tight")
    return [str(p) for p in base_dir.glob("*.png")]
