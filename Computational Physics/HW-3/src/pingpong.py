import numpy as np
from typing import Tuple, List

from .ode_solver import ODESolver

class PingPongSystem:
    """
    乒乓球受拍打的动力学系统。
    负责管理系统的常微分方程、球拍运动状态以及碰撞事件的精确检测与处理。
    """
    
    def __init__(self, g: float = 10.0, gamma: float = 0.02, A: float = 0.02, omega: float = 12.5663706):
        """
        初始化系统参数。
        
        :param g: 重力加速度
        :param gamma: 空气阻尼系数 (正比于速度)
        :param A: 球拍振动振幅
        :param omega: 球拍振动角频率
        """
        self.g = g
        self.gamma = gamma
        self.A = A
        self.omega = omega

    def dynamics(self, t: float, u: np.ndarray) -> np.ndarray:
        """
        乒乓球在空中的自然下落段的动力学方程 f(t, u)。
        状态向量 u = [y, v]^T
        
        :param t: 当前时间
        :param u: 当前状态 
        :return: 导数数组 [v, a]
        """
        y, v = u
        a = -self.g - self.gamma * v
        return np.array([v, a], dtype=float)
        
    def racket_pos(self, t: float) -> float:
        """球拍在时间 t 的位置 h(t) = A * sin(omega * t)"""
        return self.A * np.sin(self.omega * t)
        
    def racket_vel(self, t: float) -> float:
        """球拍在时间 t 的速度 v_h(t) = A * omega * cos(omega * t)"""
        return self.A * self.omega * np.cos(self.omega * t)

    def solve_with_collisions(self, t_span: Tuple[float, float], u0: np.ndarray, dt: float = 0.01) -> Tuple[np.ndarray, np.ndarray]:
        """
        求解包含突发碰撞事件的运动方程。
        碰撞判断：当 y_next <= h_next 且 相对速度 v < v_h 时发生弹性碰撞。
        使用二分法(Bisection)精确定位碰撞时间 t_c 并处理状态突变。
        
        :param t_span: 演化时间段 (t_start, t_end)
        :param u0: 初始状态 [y0, v0]
        :param dt: 标称步长
        :return: 时间数组 t_eval, 状态数组 u_eval
        """
        t_start, t_end = t_span
        n_steps = int(np.ceil((t_end - t_start) / dt))
        dt = (t_end - t_start) / n_steps
        
        t_list = [t_start]
        u_list = [u0]
        
        t = t_start
        u = np.array(u0, dtype=float)
        
        while t < t_end:
            u_next = ODESolver.rk4_step(self.dynamics, t, u, dt)
            t_next = t + dt
            
            # 使用下一拍的精确位置检测是否越界穿界
            h_next = self.racket_pos(t_next)
            
            # 碰撞条件:
            # 1. 乒乓球进入或穿过球拍面 (u_next[0] <= h_next)
            # 2. 乒乓球相对于球拍具有正在接近的趋势，以免在此处粘连
            if u_next[0] <= h_next and u[1] <= self.racket_vel(t):
                # 利用二分法精确查找碰撞发生的瞬时点 t_c
                t_L, t_R = t, t_next
                u_L = u.copy()
                
                # 设置搜索精度
                for _ in range(20):
                    t_M = (t_L + t_R) / 2.0
                    dt_M = t_M - t_L
                    u_M = ODESolver.rk4_step(self.dynamics, t_L, u_L, dt_M)
                    
                    if u_M[0] <= self.racket_pos(t_M):
                        # 碰撞在此刻之前已经发生
                        t_R = t_M
                    else:
                        # 碰撞在此刻之后才发生
                        t_L = t_M
                        u_L = u_M
                
                # 记录碰撞时刻的状态 t_c
                t_c = t_L
                u_c = u_L
                
                # 进行完全弹性碰撞更新状态
                v_h_c = self.racket_vel(t_c)
                u_c[0] = self.racket_pos(t_c) # 将球严格拉回到拍面上
                # v_after - v_h = -(v_before - v_h)
                u_c[1] = 2.0 * v_h_c - u_c[1]
                
                # 推进剩余的时间步 (从 t_c 到 t_next)
                dt_rem = t_next - t_c
                u_next = ODESolver.rk4_step(self.dynamics, t_c, u_c, dt_rem)

                # 将碰撞点也可选地加入列表，使折线更锐利清晰
                t_list.append(t_c)
                u_list.append(u_c.copy())
            
            t = t_next
            u = u_next
            
            t_list.append(t)
            u_list.append(u.copy())
            
        return np.array(t_list), np.array(u_list)
