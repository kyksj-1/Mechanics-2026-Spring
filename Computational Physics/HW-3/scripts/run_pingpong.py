import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.append(str(ROOT_DIR))

from config.config import CONFIG_DIR, OUTPUTS_DIR
from src.pingpong import PingPongSystem

def run_single_simulation(config_params, y0, t_span, dt):
    """
    运行单个初值条件下的仿真并绘制结果
    """
    print(f"Running Ping-Pong simulation for y0 = {y0}...")
    
    system = PingPongSystem(
        g=config_params["g"],
        gamma=config_params["gamma"],
        A=config_params["A"],
        omega=config_params["omega"]
    )
    
    u0 = [y0, 0.0]
    t_eval, u_eval = system.solve_with_collisions(t_span, u0, dt=dt)
    
    y_eval = u_eval[:, 0]
    h_eval = system.racket_pos(t_eval)
    
    # 画图：t in (0, 10) 以及 t in (990, 1000)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # 选择 t in [0, 10]
    mask1 = (t_eval >= 0) & (t_eval <= 10)
    ax1.plot(t_eval[mask1], y_eval[mask1], label="Ping-Pong $y(t)$", color="blue")
    ax1.plot(t_eval[mask1], h_eval[mask1], label="Racket $h(t)$", color="red", linestyle="--")
    ax1.set_title(f"Early Stage $t \in [0, 10], y_0={y0}$")
    ax1.set_xlabel("Time $t$")
    ax1.set_ylabel("Height $y$")
    ax1.grid(True)
    ax1.legend()
    
    # 选择 t in [990, 1000]
    mask2 = (t_eval >= 990) & (t_eval <= 1000)
    ax2.plot(t_eval[mask2], y_eval[mask2], label="Ping-Pong $y(t)$", color="blue")
    ax2.plot(t_eval[mask2], h_eval[mask2], label="Racket $h(t)$", color="red", linestyle="--")
    ax2.set_title(f"Late Steady Stage $t \in [990, 1000], y_0={y0}$")
    ax2.set_xlabel("Time $t$")
    ax2.set_ylabel("Height $y$")
    ax2.grid(True)
    ax2.legend()
    
    plt.tight_layout()
    output_path = OUTPUTS_DIR / f"pingpong_y0_{y0}_timeline.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")

def run_mode_analysis_simulation(config_params, y0_list, t_span, dt):
    """
    探究不同初值 y0 下的稳态运动模式
    """
    print(f"Running Ping-Pong mode analysis for different y0...")
    
    system = PingPongSystem(
        g=config_params["g"],
        gamma=config_params["gamma"],
        A=config_params["A"],
        omega=config_params["omega"]
    )
    
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ['blue', 'green', 'magenta', 'orange', 'cyan']
    
    for i, y0 in enumerate(y0_list):
        u0 = [y0, 0.0]
        # 使用较大的步长提升探究速度，由于有Bisection处理碰撞，能保持较好精度
        t_eval, u_eval = system.solve_with_collisions(t_span, u0, dt=dt)
        y_eval = u_eval[:, 0]
        
        # 只截取最后 10秒 (稳态区域) 来画图
        mask = (t_eval >= t_span[1] - 10) & (t_eval <= t_span[1])
        t_plot = t_eval[mask]
        y_plot = y_eval[mask]
        
        ax.plot(t_plot, y_plot, label=f"$y_0={y0}$", color=colors[i % len(colors)], alpha=0.8)

    # 画出球拍的高频作为对比参考
    h_plot = system.racket_pos(t_plot)
    ax.plot(t_plot, h_plot, label="Racket $h(t)$", color="red", linestyle="--", linewidth=1.5)
    
    ax.set_title(f"Modes Comparison at Steady Stage $t \in [990, 1000]$")
    ax.set_xlabel("Time $t$")
    ax.set_ylabel("Height $y$")
    ax.grid(True)
    ax.legend(loc="upper right")
    
    plt.tight_layout()
    output_path = OUTPUTS_DIR / "pingpong_modes_comparison.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")

def main():
    with open(CONFIG_DIR / "config.yaml", "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)["pingpong"]
        
    params = config["parameters"]
    y0_target = config["initial_conditions"]["y_0"]
    y0_list = config["initial_conditions"]["y_0_list"]
    
    t_span = (config["time"]["t_start"], config["time"]["t_end"])
    dt = config["time"]["dt"]
    
    # 1. 对问题 B.3.i 单独做 y0=0.3 时的双屏演化对比图
    run_single_simulation(params, y0_target, t_span, dt)
    
    # 2. 对问题 B.3.ii 做多重 y0 释放的相空间收敛模式对比
    run_mode_analysis_simulation(params, y0_list, t_span, dt)

if __name__ == "__main__":
    main()
