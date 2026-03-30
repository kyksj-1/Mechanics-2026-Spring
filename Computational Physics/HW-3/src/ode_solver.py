import numpy as np
from typing import Callable, Tuple, Optional

class ODESolver:
    """
    通用常微分方程求解器类。
    提供常见的数值积分方法，以面向对象和函数接口解耦设计。
    """
    
    @staticmethod
    def rk4_step(f: Callable, t: float, y: np.ndarray, dt: float, *args) -> np.ndarray:
        """
        单步 Runge-Kutta 4阶 (RK4) 求解。
        
        :param f: 使得 dy/dt = f(t, y, *args) 的函数
        :param t: 当前时间
        :param y: 当前状态向量
        :param dt: 时间步长
        :param args: f()所需要的额外参数
        :return: 下一步的状态向量 y_{n+1}
        """
        k1 = f(t, y, *args)
        k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1, *args)
        k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2, *args)
        k4 = f(t + dt, y + dt * k3, *args)
        
        return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    @staticmethod
    def solve_ivp_rk4(
        f: Callable, 
        t_span: Tuple[float, float], 
        y0: np.ndarray, 
        dt: float = 0.01, 
        args: Optional[tuple] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        利用经典的4阶Runge-Kutta方法求解初值问题(IVP)。
        它的接口设计受启于 scipy.integrate.solve_ivp。
        
        :param f: 右端函数, f(t, y, *args)
        :param t_span: 积分区间 (t_start, t_end)
        :param y0: 初始状态的一维数组
        :param dt: 恒定步长
        :param args: 传递给f()的可变参数元组
        :return: t_eval(时间点数组), y_eval(对应的状态数组, 形状为 (len(t_eval), len(y0)))
        """
        if args is None:
            args = ()
            
        t_start, t_end = t_span
        # 根据步长生成等距时间网格
        n_steps = int(np.ceil((t_end - t_start) / dt))
        t_eval = np.linspace(t_start, t_end, n_steps + 1)
        
        # 初始化解数组
        y_eval = np.zeros((n_steps + 1, len(y0)))
        y_eval[0] = y0
        
        # 实际积分步长（匹配linspace的步长，避免微小误差）
        actual_dt = t_eval[1] - t_eval[0]
        
        # 进行逐步积分
        y_current = np.array(y0, dtype=float)
        for i in range(n_steps):
            y_next = ODESolver.rk4_step(f, t_eval[i], y_current, actual_dt, *args)
            y_eval[i+1] = y_next
            y_current = y_next
            
        return t_eval, y_eval

# 也可以提供一个顶层函数简化调用
def solve_ivp(f, t_span, y0, method="RK4", dt=0.01, args=()) -> Tuple[np.ndarray, np.ndarray]:
    """快捷入口"""
    if method == "RK4":
        return ODESolver.solve_ivp_rk4(f, t_span, y0, dt, args)
    else:
        raise NotImplementedError(f"Method '{method}' is not implemented yet.")
