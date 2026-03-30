import numpy as np
from typing import Tuple

class KapitzaPendulum:
    """
    Kapitza 摆物理系统定义。
    提供可以直接给 ODE 求解器使用的含时微分方程右端项 f(t, u, p)。
    """
    
    def __init__(self, l: float = 1.0, m: float = 1.0, g: float = 1.0, a: float = 0.1, omega: float = 5.0):
        """
        初始化系统参数。
        
        :param l: 摆长
        :param m: 小球质量
        :param g: 重力加速度
        :param a: 驱动底座振幅
        :param omega: 驱动圆频率
        """
        self.l = l
        self.m = m
        self.g = g
        self.a = a
        self.omega = omega

    def dynamics(self, t: float, u: np.ndarray) -> np.ndarray:
        """
        动力学方程 f(t, u, p)。
        系统状态 u = [theta, dot_theta]^T
        
        :param t: 当前时刻
        :param u: 状态向量数组 [theta, dot_theta]
        :return: 导数数组 [dot_theta, ddot_theta]
        """
        theta = u[0]
        dot_theta = u[1]
        
        # ddot_theta = - (g/l - (a*omega^2/l) * cos(omega*t)) * sin(theta)
        term1 = self.g / self.l
        term2 = (self.a * self.omega**2 / self.l) * np.cos(self.omega * t)
        
        ddot_theta = - (term1 - term2) * np.sin(theta)
        
        return np.array([dot_theta, ddot_theta], dtype=float)

    def get_cartesian_coords(self, t_eval: np.ndarray, y_eval: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        将极坐标转换为笛卡尔坐标，方便可视化轨迹或做成动画。
        注意：定义 theta 偏离 y 轴负方向，所以球在 (0, -l) 未必静止(由于外加振荡)。
        
        :param t_eval: 时间数组
        :param y_eval: 状态流 ([theta, dot_theta] 的数组)
        :return: x, y 的轨迹数组
        """
        theta = y_eval[:, 0]
        # x = l * sin(theta)
        x = self.l * np.sin(theta)
        # y_0(t) = a * cos(omega * t) 作为悬挂点， y = y_0(t) - l * cos(theta)
        y = self.a * np.cos(self.omega * t_eval) - self.l * np.cos(theta)
        
        return x, y
