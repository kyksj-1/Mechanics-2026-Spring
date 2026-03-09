# scripts/question_1_run.py
import sys
import os
# 确保可以导入 src 目录下的模块
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from src.basic_component import GridPlaygroundContainer

def run_density_evolution(L: int = 32, t_steps: int = 20000, save: bool = True, show: bool = True):
    # 1. 实例化核心引擎与沙堆网格
    container = GridPlaygroundContainer(L)
    playground = np.zeros((L, L), dtype=int)
    densities = np.zeros(t_steps)
    
    print(f"开始演化: 尺寸 L={L}, 总步数={t_steps}")
    
    # 2. 演化循环
    for t in range(t_steps):
        playground, _, density = container.one_step(playground)
        densities[t] = density
        
        if (t + 1) % 5000 == 0:
            print(f"已完成 {t + 1} 步, 当前密度: {density:.4f}")
            
    # 3. 绘图渲染
    plt.figure(figsize=(10, 6))
    plt.plot(range(t_steps), densities, color='#1f77b4', alpha=0.9, linewidth=1.5, label=f'L={L} Simulation')
    
    # 二维BTW模型的理论稳态临界密度近似为 2.125
    plt.axhline(y=2.125, color='red', linestyle='--', linewidth=2, label='Theoretical Critical Density $n_c \\approx 2.125$')
    
    plt.title("Q2.1: Evolution of Average Density $n(t)$", fontsize=14, fontweight='bold')
    plt.xlabel("Time steps $t$", fontsize=12)
    plt.ylabel("Average density $n(t)$", fontsize=12)
    plt.legend(fontsize=12, loc='lower right')
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()

    if save:
        root = Path(__file__).resolve().parents[1]
        plots_dir = root / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        out_path = plots_dir / f"q2_1_density_L{L}_t{t_steps}.png"
        plt.savefig(out_path, dpi=200, bbox_inches="tight")
        print(f"图片已保存到: {out_path}")

    if show:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--L", type=int, default=32)
    parser.add_argument("--t_steps", type=int, default=20000)
    parser.add_argument("--no_save", action="store_true")
    parser.add_argument("--no_show", action="store_true")
    args = parser.parse_args()
    run_density_evolution(L=args.L, t_steps=args.t_steps, save=(not args.no_save), show=(not args.no_show))
