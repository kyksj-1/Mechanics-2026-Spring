# Computational Physics - Homework 3

## 报告信息
- **执行人**：kyksj-1
- **环境**：`conda` (`research_env`)
- **运行时间**：2026年3月
- **任务目标**：常微分方程数值求解（Runge-Kutta 法应用）与混沌/动力学系统分析

---

## 问题A：Kapitza 摆 (Kapitza Pendulum)

### A.1 系统建模与运动方程

取 $\theta$ 作为广义坐标（摆偏离 $y$ 轴负方向的角度），向下为 $\theta=0$，向上 $\theta=\pi$。
摆锤的坐标为：
$$
x = l \sin \theta \\
y = a \cos(\omega t) - l \cos \theta
$$
对其求导，得到速度的平方：
$$
\dot{x}^2 + \dot{y}^2 = l^2 \dot{\theta}^2 - 2 a l \omega \sin(\omega t) \sin\theta \, \dot{\theta} + a^2 \omega^2 \sin^2(\omega t)
$$

因此，系统的拉格朗日量 $L = T - V$ 为：
$$
L = \frac{1}{2} m \left[ l^2 \dot{\theta}^2 - 2 a l \omega \sin(\omega t) \sin\theta \, \dot{\theta} + a^2 \omega^2 \sin^2(\omega t) \right] - mg(a \cos\omega t - l \cos\theta)
$$

应用欧拉-拉格朗日方程 $\frac{d}{dt} \left( \frac{\partial L}{\partial \dot{\theta}} \right) - \frac{\partial L}{\partial \theta} = 0$ 进行推导：
1. $\frac{\partial L}{\partial \dot{\theta}} = m l^2 \dot{\theta} - m a l \omega \sin(\omega t) \sin\theta$
2. $\frac{d}{dt} \left( \frac{\partial L}{\partial \dot{\theta}} \right) = m l^2 \ddot{\theta} - m a l \omega^2 \cos(\omega t) \sin\theta - m a l \omega \sin(\omega t) \cos\theta \dot{\theta}$
3. $\frac{\partial L}{\partial \theta} = - m a l \omega \sin(\omega t) \cos\theta \dot{\theta} - mgl \sin \theta$

两者相减后化简得到系统的运动方程：
$$
m l^2 \ddot{\theta} - m a l \omega^2 \cos(\omega t) \sin\theta + mgl \sin \theta = 0
$$
化简得到角加速度：
$$
\ddot{\theta} = - \left( \frac{g}{l} - \frac{a \omega^2}{l} \cos(\omega t) \right) \sin \theta
$$

设系统的广义状态变量为 $u = [\theta, \dot{\theta}]^T$，以及系统参数 $p = \{l, m, g, a, \omega\}$，则运动方程可以显式地转换为如下的一阶偏微分方程组的向量形式：
$$
\frac{d}{dt}u(t) = f(u,t,p) = 
\begin{bmatrix}
u_1 \\
- \left( \frac{g}{l} - \frac{a \omega^2}{l} \cos(\omega t) \right) \sin u_0
\end{bmatrix}
$$

### A.2 Runge-Kutta 等数值方法实现

在 `src/ode_solver.py` 中，我们编写了一个经典的 4阶 Runge-Kutta 求解器（RK4），接口受 SciPy 的 `solve_ivp` 设计启发，使得任何满足 $f(t, y, *args)$ 的一阶ODE组在此接口中都可以直接调用。求解器的核心代码摘录如下：
```python
def rk4_step(f, t, y, dt, *args):
    k1 = f(t, y, *args)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1, *args)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2, *args)
    k4 = f(t + dt, y + dt * k3, *args)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
```
> **显式写出动力学 $f(u,t,p)$ 的意义：** 任何高阶常微分方程或方程组都可以并应当转化为第一阶常微分方程的向量流形式（如欧拉观点）。因为标准的数值方法（如RK4）依赖于根据当前状态 $(t, u)$ 前向计算局部切向量向来逼近下一步轨迹，我们只需知道系统当前流形的切向量。这种抽象让同样的 ODE 求解代码具备了极高的**通用性**。

### A.3 物理现象发现

针对 $l=m=g=1, a=0.1$ 和初始条件 $\theta(0) = \frac{4}{5}\pi, \dot{\theta}(0)=0$：
我们使用编写好的求解器代入不同频率 $\omega$ 进行积分（见 \`outputs/kapitza_omega_comparison.png\`），可以观察到截然不同的动力学现象：
1. **$\omega = 5$**：运动极不稳定，摆锤不仅会迅速掉落，还会受强迫振动发生**绕轴翻旋**（$\theta$ 随时间向负方向不断增加超过 $2\pi$ 并持续震荡漂移）。
2. **$\omega = 10$**：高频分量不足以维持小球于倒置态，摆锤向下坠落，并越过最低点在 $\theta \approx 0$（即自然悬垂位置）附近进行长期的**周期振荡**（围绕下稳定点）。
3. **$\omega = 20$**：一个奇妙的现象发生了。尽管摆的初始位置 $\frac{4}{5}\pi$ 受重力会往下掉，但极高频的底座振荡导致摆锤没能掉下，而是反而向上靠拢，然后在 $\theta \approx \pi$（即**倒置/倒立状态**）附近稳定地振荡！

### A.4 理论合理解释

上述现象可以用**有效势能法 (Effective Potential)** 进行理论解释：
Kapitza摆存在快变量（底座的高频震动产生）和慢变量（宏观摆动）。在 $\omega \gg \sqrt{g/l}$ 且振幅微小的条件下，将角度分离为慢速漂移分量 $\Theta$ 与高频微小颤动分量 $\xi$，经过平均法积分掉高频项后，系统有效势能 $V_{\text{eff}}(\theta)$ 近似由重力势能和高频动能的皮动势加成构成：
$$
V_{\text{eff}}(\theta) \approx mgl (1 - \cos\theta) + \frac{m a^2 \omega^2}{4} \sin^2\theta
$$

我们需要分析 $\theta = \pi$ 处的稳定性（在最高点 $\cos\pi = -1, \sin\pi = 0$）：
其对 $\theta$ 的二阶导数为：
$$
\left. \frac{d^2 V_{\text{eff}}}{d\theta^2} \right|_{\theta=\pi} = mgl\cos\pi + \frac{m a^2 \omega^2}{2} (\cos^2\pi - \sin^2\pi) = -mgl + \frac{m a^2 \omega^2}{2}
$$
当二阶导数 $>0$ 时，这个倒置平衡点就会从“不稳定鞍点”变为“稳定极小值点”。即：
$$
\frac{a^2 \omega^2}{2} > gl \implies \omega > \sqrt{\frac{2gl}{a^2}} = \frac{\sqrt{2} \times 1}{0.1} \approx 14.14 \text{ rad/s}
$$

*当 $\omega = 5, 10$ 时，$\omega < 14.14$，此条件未满足，有效势无法形成倒置点的势阱，重力主导摆锤下落；*
*当 $\omega = 20$ 时，$\omega > 14.14$，倒置位置有效势能存在一个局域稳定势阱，动能带来的高频压制力超过了重力倾覆力矩，从而把小摆锤稳稳地支撑在了空中。这就解释了为什么 $\omega=20$ 下摆会悬停并振荡。*
