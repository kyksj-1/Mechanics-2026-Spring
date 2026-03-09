# scripts/question_3_4_run.py
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import argparse
import matplotlib.pyplot as plt
from pathlib import Path
from src.container_parallel import ContainerParallel, log_binning, fit_powerlaw

def run_scaling_analysis(
    Ls: list[int] = [16, 32, 64, 128],
    t_steps: int = 150000,
    binning_a: float = 1.2,
    tau: float = None,
    D: float = 2.0,
    save: bool = True,
    show: bool = True,
):
    # 1. 多进程并发计算
    manager = ContainerParallel(Ls, t_steps)
    manager.run_parallel()

    # 若 tau 未指定，从最大 L 的数据自动拟合
    if tau is None:
        max_L = max(Ls)
        s_fit, p_fit = log_binning(manager.results[max_L], a=binning_a)
        tau, _, _ = fit_powerlaw(s_fit, p_fit)
        print(f"从 L={max_L} 数据自动拟合得到 tau = {tau:.3f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for idx, L in enumerate(Ls):
        scores = manager.results[L]
        s_centers, prob = log_binning(scores, a=binning_a)
        
        # --- 图 1: 有限尺寸效应 ---
        ax1.loglog(s_centers, prob, marker='o', linestyle='-', color=colors[idx],
                   markersize=4, alpha=0.8, label=f'L={L}')
        
        # --- 图 2: 标度折叠 (Data Collapse) ---
        # x轴无量纲化： s / L^D
        x_scaled = s_centers / (L ** D)
        # y轴重标度： P(s) * s^tau
        y_scaled = prob * (s_centers ** tau)
        
        ax2.loglog(x_scaled, y_scaled, marker='o', linestyle='-', color=colors[idx],
                   markersize=4, alpha=0.8, label=f'L={L}')
                   
    # 图 1 装饰
    ax1.set_title("Q3: Finite-Size Effect on $P(s)$", fontsize=14, fontweight='bold')
    ax1.set_xlabel("Avalanche Score $s$", fontsize=12)
    ax1.set_ylabel("Probability Density $P(s)$", fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, which="both", ls=":", alpha=0.6)
    
    # 图 2 装饰
    ax2.set_title("Q4: Data Collapse & Universal Scaling", fontsize=14, fontweight='bold')
    ax2.set_xlabel("Scaled Size $s / L^D$", fontsize=12)
    ax2.set_ylabel("Scaled Probability $P(s) s^{\\tau}$", fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, which="both", ls=":", alpha=0.6)
    
    plt.tight_layout()
    if save:
        root = Path(__file__).resolve().parents[1]
        plots_dir = root / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        out_path = plots_dir / f"q2_3_4_scaling_Ls{'-'.join(map(str, Ls))}_t{t_steps}.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        print(f"图片已保存到: {out_path}")

    if show:
        plt.show()

# Windows/macOS 下多进程必须在 __main__ 内触发
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ls", type=str, default="16,32,64,128")
    parser.add_argument("--t_steps", type=int, default=150000)
    parser.add_argument("--binning_a", type=float, default=1.2)
    parser.add_argument("--tau", type=float, default=None)
    parser.add_argument("--D", type=float, default=2.0)
    parser.add_argument("--no_save", action="store_true")
    parser.add_argument("--no_show", action="store_true")
    args = parser.parse_args()

    Ls = [int(x.strip()) for x in args.Ls.split(",") if x.strip()]
    run_scaling_analysis(
        Ls=Ls,
        t_steps=args.t_steps,
        binning_a=args.binning_a,
        tau=args.tau,
        D=args.D,
        save=(not args.no_save),
        show=(not args.no_show),
    )
