import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# 添加 src 和 config 到模块检索路径
ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.append(str(ROOT_DIR))

from config.config import CONFIG_DIR, OUTPUTS_DIR
from src.kapitza import KapitzaPendulum
from src.ode_solver import solve_ivp

def main():
    # 1. 加载配置
    with open(CONFIG_DIR / "config.yaml", "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)["kapitza"]
        
    params = config["parameters"]
    omegas = config["omegas"]
    t_start = config["time"]["t_start"]
    t_end = config["time"]["t_end"]
    dt = config["time"]["dt"]
    theta_0 = config["initial_conditions"]["theta_0"]
    omega_0 = config["initial_conditions"]["omega_0"]
    
    y0 = np.array([theta_0, omega_0])
    t_span = (t_start, t_end)
    
    # 画图对象初始化
    plt.figure(figsize=(15, 5))
    
    # 为每一个 omega 进行求解
    for idx, w in enumerate(omegas):
        print(f"Solving Kapitza Pendulum for omega = {w} ...")
        
        # 构建物理系统
        system = KapitzaPendulum(
            l=params["l"],
            m=params["m"],
            g=params["g"],
            a=params["a"],
            omega=w
        )
        
        # 封装右端函数，ODE求解器只需要 f(t, u)
        def f_system(t, u):
            return system.dynamics(t, u)
            
        # 调用自行编写的 RK4 求解器
        t_eval, y_eval = solve_ivp(f_system, t_span, y0, method="RK4", dt=dt)
        theta_eval = y_eval[:, 0]
        
        # 结果可视化 - 画出 theta 随时间变化的图像
        plt.subplot(1, 3, idx+1)
        plt.plot(t_eval, theta_eval, label=rf"$\omega={w}$", color="blue")
        
        # 画两条常见的参考线: theta=0(最低点) 和 theta=pi(最高点)
        plt.axhline(0, color="green", linestyle="--", alpha=0.5, label=r"Lower eq ($\theta=0$)")
        plt.axhline(np.pi, color="red", linestyle="--", alpha=0.5, label=r"Inverted eq ($\theta=\pi$)")
        
        plt.title(rf"Kapitza Pendulum: $\omega={w}$")
        plt.xlabel("Time $t$")
        plt.ylabel(r"Angle $\theta(t)$ (rad)")
        plt.ylim(-0.5, np.pi + 1.0)
        plt.legend(loc="upper right")
        plt.grid(True)

    plt.tight_layout()
    output_path = OUTPUTS_DIR / "kapitza_omega_comparison.png"
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved successfully at {output_path}")

if __name__ == "__main__":
    main()
