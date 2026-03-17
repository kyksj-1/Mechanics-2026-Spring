# 第9章 应用专题（上）：量子力学、优化与深度学习

> **作者**：kyksj-1
> **风格致敬**：Gilbert Strang × 3Blue1Brown

---

## 本章导读

前八章系统建立了特征值、对角化、二次型、正交性、正定性等工具。本章要回答一个终极问题：

> **这些工具在真实世界中究竟有什么用？**

我们将穿越五大领域——量子力学、优化与输运、深度学习、金融、刚体力学——展示线性代数如何成为现代科学与工程的**通用语言**。本篇（上）聚焦前三个领域，下篇聚焦金融与刚体力学，并给出全景总结。

每个领域都给出**非入门级**的例子和完整的代码实现。

---

## 9.1 量子力学中的应用

### 9.1.1 可观测量与厄米算符

在量子力学中，物理可观测量（能量、动量、角动量等）对应**厄米算符**（Hermitian operator），即满足 $\hat{A}^\dagger = \hat{A}$ 的算符。

在有限维空间中，厄米算符就是**厄米矩阵**（$A^\dagger = A$，实矩阵退化为对称矩阵）。

**谱定理的物理意义**：

| 数学 | 物理 |
|------|------|
| 特征值 $\lambda_i$ | 可观测量的**测量值**（一定是实数） |
| 特征向量 $\|n\rangle$ | 对应测量值的**本征态** |
| 谱分解 $\hat{A} = \sum \lambda_i \|i\rangle\langle i\|$ | 可观测量在本征态上的展开 |
| 概率 $|\langle i|\psi\rangle|^2$ | 测量得到 $\lambda_i$ 的概率 |

### 9.1.2 例题：三维耦合势场中的能级

考虑一个粒子处于三维耦合谐振势中：

$$
V(x, y, z) = \frac{1}{2}(k_{xx}x^2 + k_{yy}y^2 + k_{zz}z^2 + 2k_{xy}xy + 2k_{xz}xz + 2k_{yz}yz)
$$

这个势能是一个**二次型**：$V = \frac{1}{2}\mathbf{r}^T K \mathbf{r}$，其中 $K$ 是刚度矩阵（对称正定）。

**核心洞见**：通过正交对角化 $K = Q\Lambda Q^T$，可以将耦合的三维问题分解为三个**独立的**一维谐振子。

在新坐标 $\boldsymbol{\xi} = Q^T \mathbf{r}$ 下：

$$
V = \frac{1}{2}(\lambda_1 \xi_1^2 + \lambda_2 \xi_2^2 + \lambda_3 \xi_3^2)
$$

总能量本征值为：

$$
\boxed{E_{n_1 n_2 n_3} = \hbar\omega_1(n_1 + \tfrac{1}{2}) + \hbar\omega_2(n_2 + \tfrac{1}{2}) + \hbar\omega_3(n_3 + \tfrac{1}{2})}
$$

其中 $\omega_i = \sqrt{\lambda_i / m}$，$n_i = 0, 1, 2, \ldots$

```python
import numpy as np
from itertools import product

def coupled_harmonic_oscillator_3d(K, m=1.0, hbar=1.0, n_levels=10):
    """
    求解三维耦合谐振势的能级。

    参数:
        K: 3x3 刚度矩阵（对称正定）
        m: 粒子质量
        hbar: 约化普朗克常数
        n_levels: 显示的能级数
    """
    # 正交对角化
    eigenvalues, Q = np.linalg.eigh(K)
    omega = np.sqrt(eigenvalues / m)

    print(f"刚度矩阵 K:\n{K}")
    print(f"特征值 (k_1, k_2, k_3): {eigenvalues}")
    print(f"主轴频率 (omega_1, omega_2, omega_3): {omega}")
    print(f"旋转矩阵 Q:\n{np.round(Q, 4)}")

    # 计算能级
    max_n = 5
    levels = []
    for n1, n2, n3 in product(range(max_n), repeat=3):
        E = hbar * (omega[0]*(n1+0.5) + omega[1]*(n2+0.5) + omega[2]*(n3+0.5))
        levels.append((E, n1, n2, n3))

    levels.sort()
    print(f"\n前 {n_levels} 个能级:")
    for i, (E, n1, n2, n3) in enumerate(levels[:n_levels]):
        print(f"  E = {E:.4f}, (n1,n2,n3) = ({n1},{n2},{n3})")

    return eigenvalues, Q, levels


# 具体例子：有耦合项的势场
K = np.array([[5.0, 1.0, 0.5],
              [1.0, 3.0, 0.3],
              [0.5, 0.3, 4.0]])
coupled_harmonic_oscillator_3d(K)
```

### 9.1.3 自旋-1/2 系统与 Pauli 矩阵

自旋-1/2 粒子的自旋算符由 **Pauli 矩阵**表示：

$$
\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
\sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad
\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$

每个 Pauli 矩阵的特征值都是 $\pm 1$（对应自旋向上/向下），互相之间的**对易关系**决定了哪些可观测量可以同时测量。

$$
[\sigma_x, \sigma_y] = 2i\sigma_z \neq 0
$$

$\sigma_x$ 和 $\sigma_y$ 不可交换 → 不可同时对角化 → 不可同时测量（不确定性原理的矩阵本质）。

### 9.1.4 量子态的演化与矩阵指数

Schrodinger 方程的形式解：

$$
|\psi(t)\rangle = e^{-i\hat{H}t/\hbar}|\psi(0)\rangle
$$

如果 $\hat{H}$ 可对角化：$H = U\Lambda U^\dagger$，则：

$$
e^{-iHt/\hbar} = U \, \text{diag}(e^{-i\lambda_1 t/\hbar}, \ldots, e^{-i\lambda_n t/\hbar}) \, U^\dagger
$$

这就是对角化在量子力学中的直接应用。

### 9.1.5 双自旋耦合系统：4×4 矩阵的完整求解

两个自旋-1/2 粒子通过交换相互作用耦合，Hamilton 量为：

$$
\hat{H} = J\,\boldsymbol{\sigma}_1 \cdot \boldsymbol{\sigma}_2 = J(\sigma_x^{(1)}\sigma_x^{(2)} + \sigma_y^{(1)}\sigma_y^{(2)} + \sigma_z^{(1)}\sigma_z^{(2)})
$$

其中 $J$ 是耦合常数。在 $\{|\uparrow\uparrow\rangle, |\uparrow\downarrow\rangle, |\downarrow\uparrow\rangle, |\downarrow\downarrow\rangle\}$ 基底下，用 Kronecker 积构造 $4 \times 4$ 矩阵：

$$
\sigma_x^{(1)}\sigma_x^{(2)} = \sigma_x \otimes \sigma_x = \begin{pmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \end{pmatrix}
$$

类似地计算 $\sigma_y \otimes \sigma_y$ 和 $\sigma_z \otimes \sigma_z$，求和得到：

$$
\boxed{H = J\begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & -1 & 2 & 0 \\ 0 & 2 & -1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}}
$$

**观察**：$H$ 是**块对角**的！$|\uparrow\uparrow\rangle$ 和 $|\downarrow\downarrow\rangle$ 各自是 $1 \times 1$ 块（不与其他态耦合），而 $|\uparrow\downarrow\rangle$ 和 $|\downarrow\uparrow\rangle$ 构成一个 $2 \times 2$ 块——这正是 Ch8 分块矩阵思想的体现。

对 $2 \times 2$ 块 $\begin{pmatrix} -1 & 2 \\ 2 & -1 \end{pmatrix}$ 对角化，特征值为 $-3$ 和 $1$，特征向量为 $\frac{1}{\sqrt{2}}(|\uparrow\downarrow\rangle \mp |\downarrow\uparrow\rangle)$。

**完整能谱**：

| 态 | 能量 | 自旋量子数 $S$ | 简并度 | 物理名称 |
|----|------|---------------|--------|---------|
| $\frac{1}{\sqrt{2}}(|\uparrow\downarrow\rangle - |\downarrow\uparrow\rangle)$ | $-3J$ | $S=0$ | 1 | **单态**（singlet） |
| $|\uparrow\uparrow\rangle, \frac{1}{\sqrt{2}}(|\uparrow\downarrow\rangle + |\downarrow\uparrow\rangle), |\downarrow\downarrow\rangle$ | $J$ | $S=1$ | 3 | **三重态**（triplet） |

当 $J > 0$（反铁磁耦合），单态能量最低——两个自旋倾向于反平行排列。当 $J < 0$（铁磁耦合），三重态能量最低——自旋倾向于平行排列。

> **深层联系**：单态/三重态的分裂正是量子力学中**交换对称性**的体现。$H$ 与总自旋算符 $\hat{S}^2$ 对易（可同时对角化），所以 $S$ 是好量子数。这就是 Ch2 "可同时对角化 $\Leftrightarrow$ 算符对易"的物理实例。

```python
import numpy as np

def two_spin_system(J=1.0):
    """
    求解双自旋-1/2 耦合系统。

    参数:
        J: 耦合常数（J>0 反铁磁, J<0 铁磁）
    """
    # Pauli 矩阵
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)

    # Kronecker 积构造 4x4 矩阵
    H = J * (np.kron(sx, sx) + np.kron(sy, sy) + np.kron(sz, sz))

    print(f"Hamilton 量 H (J={J}):")
    print(np.round(H.real, 3))

    # 对角化
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    print(f"\n能量本征值: {eigenvalues.real}")

    # 总自旋算符 S^2
    S = [np.kron(s, I2) + np.kron(I2, s) for s in [sx, sy, sz]]
    S2 = sum(Si @ Si for Si in S)

    print(f"\n本征态分析:")
    basis_labels = ['|↑↑>', '|↑↓>', '|↓↑>', '|↓↓>']
    for i, (E, v) in enumerate(zip(eigenvalues, eigenvectors.T)):
        s2_val = (v.conj() @ S2 @ v).real
        S_quantum = (-1 + np.sqrt(1 + 4*s2_val)) / 2  # S(S+1) = s2_val
        print(f"  E = {E.real:+.1f}J, S = {S_quantum:.0f}")
        # 显示态的成分
        components = []
        for j, label in enumerate(basis_labels):
            if abs(v[j]) > 0.01:
                coeff = v[j]
                if abs(coeff.imag) < 1e-10:
                    components.append(f"{coeff.real:+.4f}{label}")
                else:
                    components.append(f"({coeff:.4f}){label}")
        print(f"    |psi> = {' '.join(components)}")

    # 时间演化：从 |↑↓> 出发
    psi0 = np.array([0, 1, 0, 0], dtype=complex)
    times = np.linspace(0, 2*np.pi / abs(J), 200)

    print(f"\n时间演化: 从 |↑↓> 出发")
    print(f"在 t = pi/(2J) 时，系统处于 |↓↑> 态（自旋翻转！）")

    # 计算 <sigma_z^(1)>(t)
    sz1 = np.kron(sz, I2)
    sz1_expect = []
    for t in times:
        U = eigenvectors @ np.diag(np.exp(-1j * eigenvalues * t)) @ eigenvectors.T.conj()
        psi_t = U @ psi0
        sz1_expect.append((psi_t.conj() @ sz1 @ psi_t).real)

    return times, sz1_expect


times, sz1 = two_spin_system(J=1.0)
```

### 9.1.6 密度矩阵与量子纠缠

当量子系统处于**混合态**（统计混合而非纯叠加）时，我们用**密度矩阵**描述：

$$
\rho = \sum_k p_k |\psi_k\rangle\langle\psi_k|
$$

密度矩阵的性质完全由线性代数决定：

| 性质 | 数学表达 | 物理含义 |
|------|---------|---------|
| 厄米性 | $\rho^\dagger = \rho$ | 可观测量的期望值是实数 |
| 正半定性 | $\rho \geq 0$ | 概率非负 |
| 迹为 1 | $\text{tr}(\rho) = 1$ | 概率归一化 |
| $\text{tr}(\rho^2) \leq 1$ | 等号当且仅当纯态 | **纯度**（purity）|

对于双自旋系统中的 Bell 态（最大纠缠态）：

$$
|\Psi^-\rangle = \frac{1}{\sqrt{2}}(|\uparrow\downarrow\rangle - |\downarrow\uparrow\rangle)
$$

对第一个自旋取偏迹（partial trace），得到**约化密度矩阵**：

$$
\rho_1 = \text{tr}_2(|\Psi^-\rangle\langle\Psi^-|) = \frac{1}{2}I_2
$$

$\rho_1$ 是最大混合态（所有特征值相等）——这意味着单看第一个自旋，它是**完全随机**的。但整个双粒子系统是纯态！这种"部分比整体更无序"的现象就是**量子纠缠**的数学特征。

**纠缠度量——von Neumann 熵**：

$$
\boxed{S(\rho) = -\text{tr}(\rho \ln \rho) = -\sum_k \lambda_k \ln \lambda_k}
$$

其中 $\lambda_k$ 是 $\rho$ 的特征值。$S = 0$ 对应纯态，$S = \ln d$ 对应最大混合态（$d$ 维空间中所有特征值相等）。

---

## 9.2 多变量优化与输运问题

### 9.2.1 Newton 法的完整推导

考虑无约束优化 $\min_\mathbf{x} f(\mathbf{x})$。在当前点 $\mathbf{x}_k$ 处 Taylor 展开到二阶：

$$
f(\mathbf{x}_k + \delta) \approx f(\mathbf{x}_k) + \nabla f^T \delta + \frac{1}{2}\delta^T H \delta
$$

令关于 $\delta$ 的导数为零：$H\delta = -\nabla f$，得到 Newton 步长：

$$
\boxed{\delta^* = -H^{-1}\nabla f}
$$

**正定性的角色**：
- $H$ 正定 → $\delta^*$ 是下降方向 → Newton 法向极小值前进
- $H$ 不定 → $\delta^*$ 可能指向鞍点方向 → 需要修正（如 trust region）

### 9.2.2 实战：Rosenbrock 函数优化

Rosenbrock 函数是一个经典的难优化问题：

$$
f(x, y) = (1-x)^2 + 100(y - x^2)^2
$$

最小值在 $(1, 1)$，但函数有一个狭长的"香蕉形"谷底，条件数极大。

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def rosenbrock(x):
    """Rosenbrock 函数"""
    return (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

def rosenbrock_grad(x):
    """Rosenbrock 梯度"""
    dfdx = -2*(1 - x[0]) - 400*x[0]*(x[1] - x[0]**2)
    dfdy = 200*(x[1] - x[0]**2)
    return np.array([dfdx, dfdy])

def rosenbrock_hessian(x):
    """Rosenbrock Hessian"""
    H = np.array([
        [2 - 400*(x[1] - 3*x[0]**2), -400*x[0]],
        [-400*x[0], 200]
    ])
    return H

# 比较不同优化方法
x0 = np.array([-1.5, 1.5])
methods = ['Nelder-Mead', 'CG', 'BFGS', 'Newton-CG']
results = {}

for method in methods:
    if method == 'Newton-CG':
        res = minimize(rosenbrock, x0, method=method,
                      jac=rosenbrock_grad, hess=rosenbrock_hessian)
    elif method in ['CG', 'BFGS']:
        res = minimize(rosenbrock, x0, method=method, jac=rosenbrock_grad)
    else:
        res = minimize(rosenbrock, x0, method=method)
    results[method] = res
    print(f"{method:15s}: x* = [{res.x[0]:.6f}, {res.x[1]:.6f}], "
          f"f(x*) = {res.fun:.2e}, nfev = {res.nfev}")

# 可视化
fig, ax = plt.subplots(figsize=(10, 8))
x_range = np.linspace(-2, 2, 200)
y_range = np.linspace(-1, 3, 200)
X, Y = np.meshgrid(x_range, y_range)
Z = (1-X)**2 + 100*(Y-X**2)**2

ax.contour(X, Y, np.log10(Z + 1), levels=30, cmap='viridis')
ax.plot(1, 1, 'r*', markersize=20, label='Minimum (1,1)')
ax.plot(x0[0], x0[1], 'ko', markersize=10, label='Start')

# 画 Hessian 在最优点的等高线（椭圆）
H_opt = rosenbrock_hessian([1, 1])
eigenvalues_H = np.linalg.eigvalsh(H_opt)
print(f"\nHessian at (1,1): eigenvalues = {eigenvalues_H}")
print(f"Condition number = {eigenvalues_H.max()/eigenvalues_H.min():.1f}")

ax.set_xlabel('x'); ax.set_ylabel('y')
ax.set_title('Rosenbrock function (log scale contours)\n'
             f'Hessian condition number at minimum: {eigenvalues_H.max()/eigenvalues_H.min():.0f}')
ax.legend()
plt.savefig('ch9_rosenbrock.png', dpi=150, bbox_inches='tight')
plt.show()
```

### 9.2.3 输运问题：最优传输

最优传输（Optimal Transport）将一个概率分布"搬运"到另一个分布，使总成本最小。

离散形式：给定源分布 $\mathbf{a}$，目标分布 $\mathbf{b}$，成本矩阵 $C$，求传输方案 $T$：

$$
\min_{T \geq 0} \langle C, T \rangle_F \quad \text{s.t.} \quad T\mathbf{1} = \mathbf{a}, \; T^T\mathbf{1} = \mathbf{b}
$$

Sinkhorn 算法利用矩阵的**交替行列缩放**高效求解正则化版本。

```python
import numpy as np

def sinkhorn(C, a, b, reg=0.1, max_iter=100):
    """
    Sinkhorn 算法求解正则化最优传输。

    参数:
        C: n x m 成本矩阵
        a: n 维源分布（和为1）
        b: m 维目标分布（和为1）
        reg: 正则化参数（熵正则）
        max_iter: 最大迭代次数

    返回:
        T: 最优传输方案
    """
    K = np.exp(-C / reg)  # Gibbs 核
    u = np.ones_like(a)

    for _ in range(max_iter):
        v = b / (K.T @ u)
        u = a / (K @ v)

    T = np.diag(u) @ K @ np.diag(v)
    cost = np.sum(T * C)
    return T, cost


# 示例：两个一维分布的传输
n = 50
x_src = np.linspace(0, 1, n)
x_tgt = np.linspace(0, 1, n)

# 源分布：双峰高斯
a = np.exp(-(x_src - 0.3)**2 / 0.01) + np.exp(-(x_src - 0.7)**2 / 0.01)
a /= a.sum()

# 目标分布：单峰高斯
b = np.exp(-(x_tgt - 0.5)**2 / 0.02)
b /= b.sum()

# 成本矩阵（欧氏距离的平方）
C = (x_src[:, None] - x_tgt[None, :])**2

T, cost = sinkhorn(C, a, b, reg=0.005)
print(f"最优传输成本: {cost:.6f}")
print(f"传输矩阵的形状: {T.shape}")
print(f"传输矩阵的行和误差: {np.max(np.abs(T.sum(axis=1) - a)):.2e}")
print(f"传输矩阵的列和误差: {np.max(np.abs(T.sum(axis=0) - b)):.2e}")
```

### 9.2.4 条件数对梯度下降的定量影响

考虑最简单的二次优化：$f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T A\mathbf{x} - \mathbf{b}^T\mathbf{x}$，其中 $A$ 对称正定。

梯度下降的迭代：$\mathbf{x}_{k+1} = \mathbf{x}_k - \alpha \nabla f = \mathbf{x}_k - \alpha(A\mathbf{x}_k - \mathbf{b})$

在最优步长 $\alpha^* = \frac{2}{\lambda_{\max} + \lambda_{\min}}$ 下，误差的收敛率为：

$$
\boxed{\frac{\|\mathbf{x}_{k+1} - \mathbf{x}^*\|}{\|\mathbf{x}_k - \mathbf{x}^*\|} \leq \frac{\kappa - 1}{\kappa + 1}}
$$

其中 $\kappa = \lambda_{\max}/\lambda_{\min}$ 是条件数。

| 条件数 $\kappa$ | 每步误差缩减率 | 达到 $10^{-6}$ 精度所需步数 |
|:-:|:-:|:-:|
| 1 | 0 | 1（一步到位） |
| 10 | 0.818 | 69 |
| 100 | 0.980 | 691 |
| 1000 | 0.998 | 6,908 |
| $10^6$ | 0.999998 | $\sim 7 \times 10^6$ |

**直觉**：条件数大 → 等高线是扁椭圆 → 梯度方向与最优方向差很远 → 走"锯齿形"路径，收敛缓慢。

```python
import numpy as np
import matplotlib.pyplot as plt

def gradient_descent_quadratic(A, b, x0, n_steps=500, step_size=None):
    """
    对二次函数 f(x) = 0.5 x^T A x - b^T x 执行梯度下降。

    参数:
        A: 对称正定矩阵
        b: 目标向量
        x0: 初始点
        n_steps: 迭代步数
        step_size: 学习率（None 则自动选择最优步长）
    返回:
        trajectory: 迭代轨迹
        errors: 误差序列
    """
    eigenvalues = np.linalg.eigvalsh(A)
    if step_size is None:
        step_size = 2.0 / (eigenvalues.max() + eigenvalues.min())

    x_star = np.linalg.solve(A, b)
    x = x0.copy()
    trajectory = [x.copy()]
    errors = [np.linalg.norm(x - x_star)]

    for _ in range(n_steps):
        grad = A @ x - b
        x = x - step_size * grad
        trajectory.append(x.copy())
        errors.append(np.linalg.norm(x - x_star))

    return np.array(trajectory), np.array(errors)


# 比较不同条件数
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for idx, kappa in enumerate([5, 50, 500]):
    A = np.diag([1.0, kappa])
    b = np.array([1.0, 1.0])
    x0 = np.array([3.0, 3.0])

    traj, errors = gradient_descent_quadratic(A, b, x0, n_steps=200)

    ax = axes[idx]
    # 等高线
    x_star = np.linalg.solve(A, b)
    xr = np.linspace(-1, 4, 100)
    yr = np.linspace(-0.5, 4, 100)
    X, Y = np.meshgrid(xr, yr)
    Z = 0.5*(A[0,0]*X**2 + A[1,1]*Y**2) - b[0]*X - b[1]*Y
    ax.contour(X, Y, Z, levels=20, cmap='Blues', alpha=0.5)

    # 轨迹
    ax.plot(traj[:, 0], traj[:, 1], 'r.-', markersize=2, lw=0.5, alpha=0.7)
    ax.plot(x0[0], x0[1], 'ko', markersize=8)
    ax.plot(x_star[0], x_star[1], 'g*', markersize=12)
    ax.set_title(f'kappa = {kappa}\nSteps to 1e-3: {np.searchsorted(-errors[::-1], -1e-3)}')
    ax.set_xlabel('x_1')
    ax.set_ylabel('x_2')
    ax.set_aspect('equal')

plt.suptitle('Gradient Descent: Condition Number Effect', fontsize=14)
plt.tight_layout()
plt.savefig('ch9_condition_number_gd.png', dpi=150, bbox_inches='tight')
plt.show()
```

### 9.2.5 预条件器：用线性代数加速优化

**核心思想**：如果条件数大是收敛慢的根源，那就用一个**预条件矩阵** $P \approx A^{-1}$ 来"改善"条件数。

预条件梯度下降：

$$
\mathbf{x}_{k+1} = \mathbf{x}_k - \alpha P\nabla f
$$

等价于对变量做替换 $\mathbf{y} = P^{-1/2}\mathbf{x}$，新问题的 Hessian 是 $P^{1/2}AP^{1/2}$。

若 $P = A^{-1}$，新 Hessian 是单位矩阵——条件数为 1，**一步收敛**！但 $P = A^{-1}$ 的计算代价等同于直接求解……

实际中使用的预条件器是"便宜但有效"的近似：

| 预条件器 | $P$ 的形式 | 条件数改善 | 代价 |
|---------|-----------|-----------|------|
| 对角预条件 | $P = \text{diag}(A)^{-1}$ | 中等 | $O(n)$ |
| 不完全 Cholesky | $P = (\tilde{L}\tilde{L}^T)^{-1}$ | 显著 | $O(n \cdot \text{nnz})$ |
| 块对角预条件 | $P = \text{diag}(A_1, \ldots, A_k)^{-1}$ | 视块结构而定 | $O(n^3/k^2)$ |

> **与 Ch7、Ch8 的联系**：Cholesky 预条件器来自 Ch7 的正定矩阵分解，块对角预条件器来自 Ch8 的分块思想。优化领域最前沿的预条件技术，本质上就是在用线性代数的工具改善线性代数的问题。

### 9.2.6 约束优化与 KKT 条件

带等式约束的优化 $\min f(\mathbf{x})$ s.t. $\mathbf{g}(\mathbf{x}) = 0$ 的 KKT 系统是一个**鞍点系统**：

$$
\begin{pmatrix} H & J^T \\ J & 0 \end{pmatrix}\begin{pmatrix} \delta\mathbf{x} \\ \boldsymbol{\lambda} \end{pmatrix} = -\begin{pmatrix} \nabla f + J^T\boldsymbol{\lambda} \\ \mathbf{g} \end{pmatrix}
$$

其中 $J$ 是约束的 Jacobi 矩阵。这个系数矩阵是**不定**的（有正有负特征值），但具有特殊的块结构——Ch8 的 Schur 补在这里大显身手：

$$
S = -JH^{-1}J^T
$$

消去 $\delta\mathbf{x}$ 后，只需在约束空间中求解一个关于 $\boldsymbol{\lambda}$ 的小规模系统。

---

## 9.3 深度学习中的应用

### 9.3.1 注意力机制 (Attention) 的线性代数本质

Transformer 的核心——**缩放点积注意力**：

$$
\text{Attention}(Q, K, V) = \text{softmax}\left(\frac{QK^T}{\sqrt{d_k}}\right) V
$$

**线性代数视角**：
- $QK^T$：计算 Query 和 Key 之间的**内积相似度矩阵**（$n \times n$）
- softmax：将相似度转化为概率权重
- 乘以 $V$：对 Value 做**加权组合**（一种"软投影"）

注意力矩阵 $A = \text{softmax}(QK^T/\sqrt{d_k})$ 是一个**行随机矩阵**（每行和为1），它本质上是一种**自适应投影**。

```python
import numpy as np

def scaled_dot_product_attention(Q, K, V):
    """
    实现缩放点积注意力。

    参数:
        Q: (seq_len, d_k) Query 矩阵
        K: (seq_len, d_k) Key 矩阵
        V: (seq_len, d_v) Value 矩阵

    返回:
        output: (seq_len, d_v) 注意力输出
        attention_weights: (seq_len, seq_len) 注意力权重矩阵
    """
    d_k = Q.shape[-1]

    # 相似度矩阵
    scores = Q @ K.T / np.sqrt(d_k)

    # softmax（数值稳定版）
    scores_shifted = scores - scores.max(axis=-1, keepdims=True)
    exp_scores = np.exp(scores_shifted)
    attention_weights = exp_scores / exp_scores.sum(axis=-1, keepdims=True)

    # 加权组合
    output = attention_weights @ V

    return output, attention_weights


# 示例
np.random.seed(42)
seq_len, d_k, d_v = 6, 4, 4
Q = np.random.randn(seq_len, d_k)
K = np.random.randn(seq_len, d_k)
V = np.random.randn(seq_len, d_v)

output, attn = scaled_dot_product_attention(Q, K, V)
print(f"注意力权重矩阵 (每行和为1):\n{np.round(attn, 3)}")
print(f"\n行和验证: {attn.sum(axis=1)}")

# 分析注意力矩阵的谱
eigenvalues_attn = np.linalg.eigvals(attn)
print(f"\n注意力矩阵的特征值模: {np.abs(eigenvalues_attn).round(4)}")
print("（行随机矩阵的最大特征值 = 1）")
```

### 9.3.2 损失函数地形与 Hessian 分析

深度学习中损失函数 $L(\boldsymbol{\theta})$ 的地形结构由 Hessian $H = \nabla^2 L$ 决定：

| Hessian 特性 | 地形含义 | 训练影响 |
|-------------|---------|---------|
| 正定 | 局部极小值（碗底） | 稳定收敛 |
| 不定（有正有负特征值） | 鞍点 | SGD 能逃逸 |
| 接近奇异（条件数极大） | 狭长谷底 | 训练困难 |
| 大量接近零的特征值 | 平坦方向 | 多解等价 |

**现代发现**：深度网络的 Hessian 通常有大量接近零的特征值和少量大特征值，形成**"块状"谱**。这意味着损失地形在大多数方向上几乎是平的。

```python
import numpy as np

def analyze_loss_landscape(model_fn, params, loss_fn, data, n_directions=2):
    """
    沿随机方向分析损失地形（简化版）。

    参数:
        model_fn: 模型函数 model_fn(params, x) -> predictions
        params: 参数向量
        loss_fn: 损失函数 loss_fn(predictions, targets) -> scalar
        data: (x, y) 数据
        n_directions: 探测方向数
    """
    x, y = data
    n_params = len(params)

    # 生成随机方向
    directions = []
    for _ in range(n_directions):
        d = np.random.randn(n_params)
        d /= np.linalg.norm(d)
        directions.append(d)

    # 沿方向扫描
    alphas = np.linspace(-1, 1, 101)
    losses = np.zeros((n_directions, len(alphas)))

    for i, d in enumerate(directions):
        for j, alpha in enumerate(alphas):
            perturbed = params + alpha * d
            pred = model_fn(perturbed, x)
            losses[i, j] = loss_fn(pred, y)

    # 数值估计 Hessian（仅沿探测方向）
    center_loss = loss_fn(model_fn(params, x), y)
    eps = 1e-4
    hessian_diag = []
    for d in directions:
        l_plus = loss_fn(model_fn(params + eps*d, x), y)
        l_minus = loss_fn(model_fn(params - eps*d, x), y)
        h = (l_plus - 2*center_loss + l_minus) / eps**2
        hessian_diag.append(h)

    print(f"损失在当前点: {center_loss:.6f}")
    print(f"沿随机方向的曲率: {hessian_diag}")
    print("正曲率 → 碗底方向; 负曲率 → 鞍点方向; 接近零 → 平坦方向")

    return alphas, losses


# 简单示例：线性回归的损失地形
np.random.seed(42)
n, d = 50, 3
X = np.random.randn(n, d)
w_true = np.array([1.0, -2.0, 0.5])
y = X @ w_true + 0.1 * np.random.randn(n)

model = lambda w, x: x @ w
loss = lambda pred, y: np.mean((pred - y)**2)

w0 = np.zeros(d)
alphas, losses = analyze_loss_landscape(model, w0, loss, (X, y))

# 解析 Hessian
H = 2/n * X.T @ X
print(f"\n解析 Hessian 的特征值: {np.linalg.eigvalsh(H)}")
print(f"条件数: {np.linalg.cond(H):.2f}")
```

### 9.3.3 扩散模型中的噪声调度

扩散模型（如 DDPM）的核心操作是逐步加噪和去噪：

$$
\mathbf{x}_t = \sqrt{\bar{\alpha}_t}\,\mathbf{x}_0 + \sqrt{1 - \bar{\alpha}_t}\,\boldsymbol{\epsilon}
$$

协方差结构：$\text{Cov}(\mathbf{x}_t) = \bar{\alpha}_t \text{Cov}(\mathbf{x}_0) + (1-\bar{\alpha}_t)I$

这是两个正定矩阵的**凸组合**（仍然正定），其特征值随 $t$ 的变化揭示了信息损失的过程。

```python
import numpy as np
import matplotlib.pyplot as plt

def diffusion_covariance_analysis(Sigma_0, n_steps=100):
    """
    分析扩散过程中协方差矩阵特征值的演化。

    参数:
        Sigma_0: 原始数据的协方差矩阵（正定）
        n_steps: 扩散步数
    """
    # 线性噪声调度
    betas = np.linspace(1e-4, 0.02, n_steps)
    alphas = 1 - betas
    alpha_bar = np.cumprod(alphas)

    eigenvalues_0 = np.linalg.eigvalsh(Sigma_0)
    d = len(eigenvalues_0)

    # 在正交基下，Sigma_t 的特征值
    eigenvalues_t = np.zeros((n_steps, d))
    for t in range(n_steps):
        ab = alpha_bar[t]
        eigenvalues_t[t] = ab * eigenvalues_0 + (1 - ab)

    # 可视化
    plt.figure(figsize=(12, 5))

    plt.subplot(121)
    for i in range(d):
        plt.plot(eigenvalues_t[:, i], label=f'lambda_{i+1}')
    plt.axhline(y=1, color='k', ls='--', alpha=0.5, label='Pure noise (I)')
    plt.xlabel('Diffusion step t')
    plt.ylabel('Eigenvalue of Cov(x_t)')
    plt.title('Eigenvalue evolution during diffusion')
    plt.legend()

    plt.subplot(122)
    kappa = eigenvalues_t.max(axis=1) / eigenvalues_t.min(axis=1)
    plt.plot(kappa, 'r-', lw=2)
    plt.xlabel('Diffusion step t')
    plt.ylabel('Condition number')
    plt.title('Condition number -> 1 (isotropic noise)')

    plt.tight_layout()
    plt.savefig('ch9_diffusion_eigenvalues.png', dpi=150, bbox_inches='tight')
    plt.show()


# 原始数据有明显的方向性（条件数大）
Sigma_0 = np.array([[5.0, 2.0, 0.5],
                    [2.0, 2.0, 0.3],
                    [0.5, 0.3, 0.5]])
diffusion_covariance_analysis(Sigma_0)
```

### 9.3.4 权重矩阵的奇异值与梯度爆炸/消失

深度网络可以看作一连串矩阵乘法。对于 $L$ 层的线性网络（忽略激活函数）：

$$
\mathbf{y} = W_L W_{L-1} \cdots W_1 \mathbf{x}
$$

反向传播中的梯度也是矩阵连乘：

$$
\frac{\partial \mathcal{L}}{\partial \mathbf{x}} = W_1^T W_2^T \cdots W_L^T \frac{\partial \mathcal{L}}{\partial \mathbf{y}}
$$

设每层权重矩阵的最大奇异值为 $\sigma_{\max}(W_k)$，则：

$$
\left\|\frac{\partial \mathcal{L}}{\partial \mathbf{x}}\right\| \leq \prod_{k=1}^{L} \sigma_{\max}(W_k) \cdot \left\|\frac{\partial \mathcal{L}}{\partial \mathbf{y}}\right\|
$$

- 若所有 $\sigma_{\max} > 1$：梯度随层数**指数增长**（梯度爆炸）
- 若所有 $\sigma_{\max} < 1$：梯度随层数**指数衰减**（梯度消失）
- 理想情况：$\sigma_{\max} \approx 1$——这就是为什么**正交初始化**（$W$ 初始为正交矩阵）效果好

$$
\boxed{\text{梯度稳定} \iff \sigma_{\max}(W_k) \approx 1, \;\forall k}
$$

**谱归一化**（Spectral Normalization）正是基于这个观察：将每层权重除以其最大奇异值，强制 $\sigma_{\max} = 1$。这在 GAN 训练中广泛使用。

```python
import numpy as np
import matplotlib.pyplot as plt

def analyze_gradient_flow(weight_matrices, input_dim=100, n_trials=50):
    """
    分析多层网络中梯度的传播行为。

    参数:
        weight_matrices: 权重矩阵列表 [W_1, W_2, ..., W_L]
        input_dim: 输入维度
        n_trials: 随机试验次数
    """
    L = len(weight_matrices)
    grad_norms = np.zeros((n_trials, L + 1))

    for trial in range(n_trials):
        # 模拟反向传播
        g = np.random.randn(weight_matrices[-1].shape[0])
        g /= np.linalg.norm(g)
        grad_norms[trial, L] = np.linalg.norm(g)

        for k in range(L - 1, -1, -1):
            g = weight_matrices[k].T @ g
            grad_norms[trial, k] = np.linalg.norm(g)

    # 计算每层的最大奇异值
    singular_values = [np.linalg.svd(W, compute_uv=False)[0]
                       for W in weight_matrices]

    return grad_norms, singular_values


# 比较三种初始化策略
d = 100
L = 20

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
strategies = {
    'Too large (sigma_max > 1)': lambda: np.random.randn(d, d) * 1.5 / np.sqrt(d),
    'Orthogonal (sigma_max = 1)': lambda: np.linalg.qr(np.random.randn(d, d))[0],
    'Too small (sigma_max < 1)': lambda: np.random.randn(d, d) * 0.5 / np.sqrt(d),
}

for idx, (name, init_fn) in enumerate(strategies.items()):
    Ws = [init_fn() for _ in range(L)]
    grads, svs = analyze_gradient_flow(Ws)

    ax = axes[idx]
    # 梯度范数的中位数和分位数
    median = np.median(grads, axis=0)
    q25, q75 = np.percentile(grads, [25, 75], axis=0)
    layers = np.arange(L + 1)

    ax.semilogy(L - layers, median, 'b-', lw=2)
    ax.fill_between(L - layers, q25, q75, alpha=0.3)
    ax.set_xlabel('Layer (from output to input)')
    ax.set_ylabel('Gradient norm (log scale)')
    ax.set_title(f'{name}\nmax sigma_1 = {max(svs):.2f}')
    ax.set_ylim(1e-10, 1e10)

plt.suptitle(f'Gradient Flow in {L}-Layer Network (d={d})', fontsize=14)
plt.tight_layout()
plt.savefig('ch9_gradient_flow.png', dpi=150, bbox_inches='tight')
plt.show()
```

### 9.3.5 Batch Normalization 的线性代数本质

**Batch Normalization** 的核心操作是对每个 mini-batch 的激活值进行**白化**（whitening）——即使其均值为零、协方差为单位矩阵。

简化版本（忽略可学习参数 $\gamma, \beta$）：

$$
\hat{\mathbf{x}} = \Sigma^{-1/2}(\mathbf{x} - \boldsymbol{\mu})
$$

其中 $\boldsymbol{\mu}$ 和 $\Sigma$ 是 mini-batch 的均值和协方差。

**为什么这能加速训练？**

1. **条件数改善**：白化后的 Hessian 条件数趋向 1（回顾 9.2.5 预条件器的思想）
2. **去除特征间相关性**：$\Sigma^{-1/2}$ 就是对角化 + 缩放，使各方向的曲率一致
3. **Internal Covariate Shift**：每层的输入分布不再剧烈变化

实际的 BatchNorm 只做**逐通道**归一化（对角元素），而不是完整的白化（需要求 $\Sigma^{-1/2}$），因为完整白化的计算代价太高。但其思想根源是线性代数中的**协方差对角化**。

### 9.3.6 低秩近似与模型压缩

训练好的神经网络，其权重矩阵往往具有**近似低秩**结构——大部分奇异值接近零。

SVD 分解 $W = U\Sigma V^T$，保留前 $r$ 个奇异值：

$$
W \approx W_r = U_r \Sigma_r V_r^T
$$

将一个 $m \times n$ 的全连接层替换为两个小矩阵的乘积：
- $U_r\Sigma_r$：$m \times r$
- $V_r^T$：$r \times n$

参数量从 $mn$ 降至 $r(m+n)$。当 $r \ll \min(m,n)$ 时，压缩比显著。

**LoRA**（Low-Rank Adaptation）就是这个思想的升级版：微调大模型时，不修改原始权重 $W_0$，只学习一个低秩增量 $\Delta W = BA$（$B \in \mathbb{R}^{m \times r}, A \in \mathbb{R}^{r \times n}$），参数量从 $mn$ 降至 $r(m+n)$。

---

## 9.4 上篇 Key Takeaway

| 领域 | 核心矩阵/算符 | 关键线性代数操作 | 物理/工程意义 |
|------|-------------|----------------|-------------|
| 量子力学 | 厄米矩阵 $\hat{A}$ | 谱分解 | 测量值 = 特征值，本征态 = 特征向量 |
| 量子力学 | 密度矩阵 $\rho$ | 偏迹、特征值分解 | 纠缠度量、纯度 |
| 量子力学 | Kronecker 积 $H_1 \otimes I + I \otimes H_2$ | 块对角化 | 复合系统分解 |
| 优化 | Hessian $H = \nabla^2 f$ | 特征值分析 | 条件数 → 收敛速度 |
| 优化 | 预条件矩阵 $P$ | $P \approx H^{-1}$ | 改善条件数 → 加速优化 |
| 优化 | KKT 鞍点系统 | Schur 补 | 约束消元 |
| 深度学习 | 注意力矩阵 | 行随机矩阵、谱分析 | 自适应加权 |
| 深度学习 | 权重矩阵 $W$ | SVD, 奇异值分析 | 梯度稳定性、模型压缩 |
| 深度学习 | Hessian 谱 | 特征值分布 | 损失地形结构 |
| 深度学习 | 协方差 $\Sigma_t$ | 凸组合、对角化 | 扩散过程的信息流 |

---

## 习题

### 概念理解

**9.1** 解释为什么量子力学中可观测量必须用厄米矩阵表示。如果用非厄米矩阵，物理上会出什么问题？

**9.2** 在优化问题中，"条件数"和"收敛速度"是什么关系？用二次函数 $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x}$ 给出定量分析。

**9.3** Transformer 中的注意力矩阵是行随机矩阵。证明行随机矩阵的最大特征值为 1，对应的右特征向量为 $\mathbf{1}$。

### 计算练习

**9.4** 一个三维耦合谐振势 $V = 2x^2 + 3y^2 + 4z^2 + 2xy - 2yz$。
  - 写出刚度矩阵 $K$
  - 求主轴方向和主频率
  - 写出前5个能级（取 $m = \hbar = 1$）

**9.5** 对矩阵 $A = \begin{pmatrix} -1 & 2 \\ -3 & -1 \end{pmatrix}$ 代表的线性系统 $\dot{\mathbf{x}} = A\mathbf{x}$：
  - 判断稳定性
  - 画出相图（编程）
  - 求 Lyapunov 函数

**9.6** 对双自旋耦合 Hamilton 量 $H = J(\boldsymbol{\sigma}_1 \cdot \boldsymbol{\sigma}_2)$：
  - 用 Kronecker 积显式构造 $4 \times 4$ 矩阵
  - 求解全部本征值和本征态
  - 验证 $[H, \hat{S}^2] = 0$（$H$ 与总自旋算符对易）
  - 如果初始态为 $|\uparrow\downarrow\rangle$，求 $t = \pi/(2J)$ 时的态

### 思考题

**9.7** 扩散模型中，为什么协方差矩阵的条件数随扩散步数趋向 1？这对"去噪"过程意味着什么？

（提示：考虑 $\text{Cov}(\mathbf{x}_t) = \bar{\alpha}_t \Sigma_0 + (1-\bar{\alpha}_t)I$ 的特征值在 $\bar{\alpha}_t \to 0$ 时的行为。）

**9.8** 比较三种加速梯度下降的策略：(a) 动量法, (b) 预条件（对角预条件）, (c) Newton 法。从条件数改善的角度分析它们各自的优势和代价。

**9.9** 解释为什么 LoRA（Low-Rank Adaptation）能有效地微调大语言模型。从权重矩阵的谱结构出发，说明低秩假设的合理性。

（提示：训练好的权重矩阵通常有快速衰减的奇异值。微调只需调整少数几个主要方向。）

### 编程题

**9.10** 实现一个完整的量子态演化模拟器：
  - 构造 Heisenberg 自旋链 $H = J\sum_{i=1}^{N-1}\boldsymbol{\sigma}_i \cdot \boldsymbol{\sigma}_{i+1}$（$N = 6$）
  - 从 Neel 态 $|\uparrow\downarrow\uparrow\downarrow\uparrow\downarrow\rangle$ 出发，模拟时间演化
  - 计算并绘制每个位置的 $\langle\sigma_z^{(i)}\rangle(t)$ 随时间的演化（类似"自旋波传播"）
  - 计算第一个自旋的约化密度矩阵和 von Neumann 熵的时间演化

**9.11** 梯度流分析实验：
  - 构造一个 $L = 30$ 层、每层 $d = 50$ 的线性网络
  - 分别用 (a) 正态初始化 $W_{ij} \sim \mathcal{N}(0, 2/d)$, (b) 正交初始化, (c) 谱归一化后的正态初始化
  - 绘制每层梯度范数（从输出层到输入层），比较三种策略
  - 进一步：如果加入 ReLU 激活函数，结论如何变化？

---

> **下篇预告**：金融中的投资组合优化、刚体力学中的 Tennis Racket 定理，以及五大领域的全景联系图——详见第9章（下）。
