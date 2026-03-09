# scripts/question_2_run.py
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from src.container_parallel import worker_func, log_binning, fit_powerlaw

def run_single_distribution(
    L: int = 64,
    t_steps: int = 100000,
    burn_in_ratio: float = 0.2,
    binning_a: float = 1.2,
    save: bool = True,
    show: bool = True,
):
    print(f"启动单尺寸长时间演化: L={L}, 总步数={t_steps}")
    # 1. 调用 worker_func 获取有效得分数据
    _, scores = worker_func(L, t_steps, burn_in_ratio=burn_in_ratio)
    
    root = Path(__file__).resolve().parents[1]
    plots_dir = root / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # 2. 图1: 概率密度分布图（每个分数对应一个频率）
    unique_scores, counts = np.unique(scores, return_counts=True)
    freq = counts / len(scores)

    plt.figure(figsize=(8, 6))
    plt.bar(unique_scores, freq, width=1.0, color='#1f77b4', alpha=0.7, edgecolor='none')
    plt.title("Q2: Score Probability Distribution $P(s)$", fontsize=14, fontweight='bold')
    plt.xlabel("Score $s$", fontsize=12)
    plt.ylabel("Frequency $P(s)$", fontsize=12)
    plt.grid(True, axis='y', ls="--", alpha=0.5)
    plt.tight_layout()

    if save:
        out_path_raw = plots_dir / f"q2_2_Ps_dist_L{L}_t{t_steps}.png"
        plt.savefig(out_path_raw, dpi=200, bbox_inches="tight")
        print(f"图片已保存到: {out_path_raw}")

    # 3. 图2: 对数分箱 + 幂律拟合
    s_centers, prob = log_binning(scores, a=binning_a)
    if len(s_centers) == 0:
        raise RuntimeError("没有得到可绘制的分布数据（可能是样本太少或分箱参数不合适）")

    plt.figure(figsize=(8, 6))
    plt.loglog(s_centers, prob, marker='o', linestyle='-', color='#2ca02c',
               markersize=5, linewidth=1.5, label=f'Log-binned (L={L})')

    # 幂律拟合并绘制拟合线
    tau, fit_s, fit_p = fit_powerlaw(s_centers, prob)
    plt.loglog(fit_s, fit_p, 'r--', linewidth=2,
               label=f'Power-law fit: $\\tau = {tau:.2f}$')
    print(f"拟合得到幂律指数 tau = {tau:.3f}")

    plt.title("Q2: Avalanche Size Distribution $P(s)$ (Log-binned)", fontsize=14, fontweight='bold')
    plt.xlabel("Avalanche Score $s$", fontsize=12)
    plt.ylabel("Probability Density $P(s)$", fontsize=12)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=12)
    plt.tight_layout()

    if save:
        out_path_binned = plots_dir / f"q2_2_Ps_binned_L{L}_t{t_steps}.png"
        plt.savefig(out_path_binned, dpi=200, bbox_inches="tight")
        print(f"图片已保存到: {out_path_binned}")

    if show:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--L", type=int, default=64)
    parser.add_argument("--t_steps", type=int, default=100000)
    parser.add_argument("--burn_in_ratio", type=float, default=0.2)  # 热启动的比例
    parser.add_argument("--binning_a", type=float, default=1.2)  # 分箱参数
    parser.add_argument("--no_save", action="store_true")
    parser.add_argument("--no_show", action="store_true")
    args = parser.parse_args()
    run_single_distribution(
        L=args.L,
        t_steps=args.t_steps,
        burn_in_ratio=args.burn_in_ratio,
        binning_a=args.binning_a,
        save=(not args.no_save),
        show=(not args.no_show),
    )
