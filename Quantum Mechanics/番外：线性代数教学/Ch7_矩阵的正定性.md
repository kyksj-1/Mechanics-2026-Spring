# 第7章 矩阵的正定性 (Positive Definiteness)

> **作者**：kyksj-1
> **风格致敬**：Gilbert Strang × 3Blue1Brown

---

## 本章导读

正定矩阵是线性代数中最"好"的矩阵。正如 Gilbert Strang 所说：

> "正定性把特征值、对角化、二次型、能量、稳定性这些看似不同的概念，串成了一条线。"

本章将系统讲解正定性的多种等价判据，并通过具体例子展示正定性在稳定性分析、优化理论、以及各工程领域中的核心角色。

---

## 7.1 正定矩阵的等价判据

### 7.1.1 定义回顾

对称矩阵 $A \in \mathbb{R}^{n\times n}$ 是**正定的**（positive definite），若：

$$
\boxed{\mathbf{x}^T A \mathbf{x} > 0, \quad \forall \mathbf{x} \neq \mathbf{0}}
$$

### 7.1.2 五个等价条件

**定理**：以下条件互相等价：

| # | 条件 | 直觉 |
|---|------|------|
| 1 | $\mathbf{x}^TA\mathbf{x} > 0$，$\forall \mathbf{x} \neq 0$ | "能量"总为正 |
| 2 | 所有特征值 $\lambda_i > 0$ | 各方向均为"拉伸" |
| 3 | 所有顺序主子式 $\Delta_k > 0$ | Sylvester 准则 |
| 4 | 存在 Cholesky 分解 $A = LL^T$ | $L$ 下三角，对角元素为正 |
| 5 | 所有主元（高斯消元的 pivot）为正 | 消元过程不出问题 |

**证明要点**：

$(1) \Leftrightarrow (2)$：已在 Ch3（特征值判别法）证明。

$(2) \Rightarrow (3)$：$\Delta_k = \det(A_k)$（$A_k$ 是 $A$ 的 $k$ 阶顺序主子矩阵），而 $A_k$ 也是对称正定的（对 $\mathbf{x}$ 的前 $k$ 个分量限制，其余取零），所以 $\det(A_k) = \prod_{i} \lambda_i^{(k)} > 0$。

$(1) \Rightarrow (4)$：Cholesky 分解的存在性证明（构造性）

对 $A$ 做高斯消元 $A = LDL^T$（$L$ 下三角，$D$ 对角）。$A$ 正定保证所有 pivot $d_i > 0$。令 $\tilde{L} = L\sqrt{D}$，则 $A = \tilde{L}\tilde{L}^T$。$\blacksquare$

### 7.1.3 Cholesky 分解

**定理**：正定矩阵 $A$ 可以**唯一**分解为：

$$
\boxed{A = LL^T}
$$

其中 $L$ 是下三角矩阵，且对角元素为正。

**算法**（$2\times 2$ 示例）：

$$
\begin{pmatrix} a & b \\ b & c \end{pmatrix} = \begin{pmatrix} l_{11} & 0 \\ l_{21} & l_{22} \end{pmatrix}\begin{pmatrix} l_{11} & l_{21} \\ 0 & l_{22} \end{pmatrix}
$$

匹配元素：
- $l_{11}^2 = a \Rightarrow l_{11} = \sqrt{a}$
- $l_{21}l_{11} = b \Rightarrow l_{21} = b/\sqrt{a}$
- $l_{21}^2 + l_{22}^2 = c \Rightarrow l_{22} = \sqrt{c - b^2/a}$

正定保证 $a > 0$ 且 $c - b^2/a > 0$（即 $\Delta_2 > 0$），所以所有根号下都为正。

**例**：$A = \begin{pmatrix} 4 & 2 \\ 2 & 5 \end{pmatrix}$

$l_{11} = 2$，$l_{21} = 1$，$l_{22} = \sqrt{5 - 1} = 2$

$$
A = \begin{pmatrix} 2 & 0 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix} \quad \checkmark
$$

> **计算优势**：Cholesky 分解比一般 LU 分解快约一倍，且数值更稳定。它是正定矩阵求解线性方程组的首选方法。

---

## 7.2 正定性与能量

### 7.2.1 能量的二次型表示

在物理中，许多能量都是二次型的形式：

$$
E = \frac{1}{2}\mathbf{x}^T K \mathbf{x}
$$

其中 $K$ 是**刚度矩阵**（弹性势能）或**质量矩阵**（动能）。

$K$ 正定意味着**能量总为正**（除了平衡点 $\mathbf{x} = 0$）。这是物理系统**稳定**的基本条件。

### 7.2.2 弹簧系统示例

考虑两个质量块通过三根弹簧连接：

```
墙 ——k1—— m1 ——k2—— m2 ——k3—— 墙
```

势能：

$$
V = \frac{1}{2}k_1 x_1^2 + \frac{1}{2}k_2(x_2 - x_1)^2 + \frac{1}{2}k_3 x_2^2
$$

$$
= \frac{1}{2}\begin{pmatrix} x_1 & x_2 \end{pmatrix}\begin{pmatrix} k_1 + k_2 & -k_2 \\ -k_2 & k_2 + k_3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix}
$$

刚度矩阵 $K = \begin{pmatrix} k_1 + k_2 & -k_2 \\ -k_2 & k_2 + k_3 \end{pmatrix}$。

当所有 $k_i > 0$ 时，$K$ 正定（可用 Sylvester 准则验证）。这意味着平衡点 $(0, 0)$ 是**稳定的**。

---

## 7.3 正定性与稳定性

### 7.3.1 Lyapunov 稳定性理论

对线性系统 $\dot{\mathbf{x}} = A\mathbf{x}$，若存在正定矩阵 $P$ 使得 $A^TP + PA$ 负定（即 $A^TP + PA < 0$），则系统**渐近稳定**。

这就是 **Lyapunov 方程**：

$$
\boxed{A^TP + PA = -Q \quad (Q \text{ 正定})}
$$

**物理直觉**：$V(\mathbf{x}) = \mathbf{x}^TP\mathbf{x}$ 是一个"能量函数"（正定保证它像一个碗）。$\dot{V} = \mathbf{x}^T(A^TP + PA)\mathbf{x} < 0$ 意味着能量严格递减，系统一定趋向原点。

### 7.3.2 示例：阻尼振荡器

$$
\ddot{x} + 2\gamma\dot{x} + \omega_0^2 x = 0 \quad (\gamma > 0)
$$

写成一阶系统 $\dot{\mathbf{x}} = A\mathbf{x}$，其中 $\mathbf{x} = (x, \dot{x})^T$：

$$
A = \begin{pmatrix} 0 & 1 \\ -\omega_0^2 & -2\gamma \end{pmatrix}
$$

特征值：$\lambda = -\gamma \pm \sqrt{\gamma^2 - \omega_0^2}$。当 $\gamma > 0$ 时，所有特征值的实部为负，系统稳定。

选取 $P = \begin{pmatrix} \omega_0^2 & 0 \\ 0 & 1 \end{pmatrix}$（正定），可以验证 $A^TP + PA$ 负半定。

---

## 7.4 正定性与优化

### 7.4.1 凸函数与 Hessian 的正定性

**定理**：二次可微函数 $f: \mathbb{R}^n \to \mathbb{R}$ 是**严格凸**的，当且仅当其 Hessian 矩阵 $H$ 在所有点上正定。

严格凸函数有**唯一的全局最小值**，梯度下降**一定收敛**。

### 7.4.2 Newton 法中的正定性

Newton 法的更新步骤：

$$
\mathbf{x}_{k+1} = \mathbf{x}_k - H^{-1} \nabla f(\mathbf{x}_k)
$$

这要求 Hessian $H$ **可逆**。如果 $H$ 正定，$-H^{-1}\nabla f$ 保证是**下降方向**，算法收敛。如果 $H$ 不定（鞍点附近），Newton 法可能走向鞍点。

### 7.4.3 条件数：优化的难度指标

正定矩阵的**条件数**：

$$
\kappa(A) = \frac{\lambda_{\max}}{\lambda_{\min}}
$$

$\kappa$ 越大，优化越困难（等高线越"扁"，梯度下降越慢）。

| $\kappa$ | 等高线形状 | 优化难度 |
|----------|-----------|---------|
| $\approx 1$ | 圆形 | 容易 |
| $\approx 10$ | 轻微椭圆 | 较容易 |
| $\approx 100$ | 扁椭圆 | 困难 |
| $\gg 1000$ | 极扁长条 | 非常困难（ill-conditioned） |

---

## 7.5 正定矩阵的运算性质

### 7.5.1 封闭性

| 运算 | 结果 | 条件 |
|------|------|------|
| $A + B$ | 正定 | $A, B$ 正定 |
| $\alpha A$ | 正定 | $A$ 正定，$\alpha > 0$ |
| $A^{-1}$ | 正定 | $A$ 正定 |
| $A^k$ | 正定 | $A$ 正定，$k \in \mathbb{Z}^+$ |
| $B^TAB$ | 正定 | $A$ 正定，$B$ 可逆 |
| $A \otimes B$ (Kronecker) | 正定 | $A, B$ 正定 |

### 7.5.2 正定矩阵的平方根

正定矩阵 $A$ 存在唯一的正定**平方根** $A^{1/2}$，使得 $A^{1/2}A^{1/2} = A$。

构造方法：$A = Q\Lambda Q^T$，则 $A^{1/2} = Q\Lambda^{1/2}Q^T$，其中 $\Lambda^{1/2} = \text{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})$。

---

## 7.6 编程实践

### 7.6.1 正定性检测与 Cholesky 分解

```python
import numpy as np
from scipy.linalg import cholesky, cho_solve, cho_factor

def analyze_positive_definiteness(A):
    """
    全面检测矩阵正定性，使用多种方法。

    参数:
        A: n x n 对称矩阵
    """
    n = A.shape[0]
    print(f"矩阵 A ({n}x{n}):")
    print(A)

    # 方法1：特征值
    eigenvalues = np.linalg.eigvalsh(A)
    print(f"\n特征值: {eigenvalues}")
    print(f"方法1（特征值）: {'正定' if all(eigenvalues > 0) else '非正定'}")

    # 方法2：顺序主子式
    minors = []
    for k in range(1, n+1):
        minors.append(np.linalg.det(A[:k, :k]))
    print(f"\n顺序主子式: {minors}")
    print(f"方法2（Sylvester）: {'正定' if all(m > 0 for m in minors) else '非正定'}")

    # 方法3：Cholesky 分解
    try:
        L = cholesky(A, lower=True)
        print(f"\n方法3（Cholesky）: 正定")
        print(f"L =\n{L}")
        print(f"验证 LL^T =\n{L @ L.T}")
    except np.linalg.LinAlgError:
        print(f"\n方法3（Cholesky）: 非正定（分解失败）")

    # 条件数
    cond = eigenvalues.max() / eigenvalues.min() if eigenvalues.min() > 0 else float('inf')
    print(f"\n条件数 kappa = {cond:.4f}")


# ============================================================
# 示例
# ============================================================

# 正定矩阵
A_pd = np.array([[4, 2, 1],
                 [2, 5, 3],
                 [1, 3, 6]])
analyze_positive_definiteness(A_pd)

print("\n" + "="*50 + "\n")

# 不定矩阵
A_indef = np.array([[1, 2],
                    [2, 1]])
analyze_positive_definiteness(A_indef)
```

### 7.6.2 条件数对梯度下降的影响

```python
import numpy as np
import matplotlib.pyplot as plt

def gradient_descent_with_condition(kappa, n_iter=50):
    """
    展示条件数对梯度下降收敛速度的影响。

    参数:
        kappa: 目标条件数
        n_iter: 迭代次数
    """
    # 构造条件数为 kappa 的正定矩阵
    Q = np.array([[np.cos(0.3), -np.sin(0.3)],
                  [np.sin(0.3),  np.cos(0.3)]])
    D = np.diag([kappa, 1.0])
    A = Q @ D @ Q.T

    b = np.array([1.0, 2.0])
    x_opt = np.linalg.solve(A, b)

    # 梯度下降：最小化 f(x) = 0.5 * x^T A x - b^T x
    x = np.array([3.0, 3.0])
    lr = 2.0 / (kappa + 1)  # 最优学习率
    trajectory = [x.copy()]
    losses = [0.5 * x @ A @ x - b @ x]

    for _ in range(n_iter):
        grad = A @ x - b
        x = x - lr * grad
        trajectory.append(x.copy())
        losses.append(0.5 * x @ A @ x - b @ x)

    return np.array(trajectory), np.array(losses), A, x_opt


fig, axes = plt.subplots(1, 3, figsize=(18, 5))
kappas = [1.5, 10, 100]

for idx, kappa in enumerate(kappas):
    traj, losses, A, x_opt = gradient_descent_with_condition(kappa)

    ax = axes[idx]
    # 等高线
    x_range = np.linspace(-1, 4, 100)
    y_range = np.linspace(-1, 4, 100)
    X, Y = np.meshgrid(x_range, y_range)
    b = A @ x_opt
    Z = 0.5 * (A[0,0]*X**2 + 2*A[0,1]*X*Y + A[1,1]*Y**2) - b[0]*X - b[1]*Y

    ax.contour(X, Y, Z, levels=30, cmap='coolwarm', alpha=0.5)
    ax.plot(traj[:, 0], traj[:, 1], 'ko-', markersize=2, linewidth=0.8)
    ax.plot(x_opt[0], x_opt[1], 'r*', markersize=15)
    ax.set_title(f'kappa = {kappa}\n({len(traj)-1} iterations)', fontsize=12)
    ax.set_aspect('equal')
    ax.set_xlim(-1, 4); ax.set_ylim(-1, 4)

plt.suptitle('Effect of Condition Number on Gradient Descent', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('ch7_condition_number.png', dpi=150, bbox_inches='tight')
plt.show()
```

### 7.6.3 Lyapunov 稳定性验证

```python
import numpy as np
from scipy.linalg import solve_lyapunov

def lyapunov_stability_check(A):
    """
    通过 Lyapunov 方程检查线性系统的稳定性。
    求解 A^T P + P A = -Q，检查 P 是否正定。

    参数:
        A: 系统矩阵
    """
    n = A.shape[0]
    Q = np.eye(n)  # 选 Q = I

    # 求解 Lyapunov 方程 A^T P + P A + Q = 0
    P = solve_lyapunov(A.T, -Q)

    print(f"系统矩阵 A:\n{A}")
    print(f"特征值: {np.linalg.eigvals(A)}")
    print(f"\nLyapunov 矩阵 P:\n{np.round(P, 4)}")

    eig_P = np.linalg.eigvalsh(P)
    print(f"P 的特征值: {np.round(eig_P, 4)}")

    if all(eig_P > 0):
        print("P 正定 => 系统渐近稳定")
    else:
        print("P 非正定 => 无法通过此方法确认稳定性")

    # 验证
    residual = A.T @ P + P @ A + Q
    print(f"验证 A^T P + P A + Q = 0 的残差: {np.linalg.norm(residual):.2e}")


# 稳定系统
A_stable = np.array([[-1, 0.5],
                     [0, -2]])
print("=== 稳定系统 ===")
lyapunov_stability_check(A_stable)

# 不稳定系统
A_unstable = np.array([[1, 0],
                       [0, -2]])
print("\n=== 不稳定系统 ===")
lyapunov_stability_check(A_unstable)
```

---

## 7.7 Key Takeaway

| 概念 | 核心要点 |
|------|---------|
| 正定 $\Leftrightarrow$ 所有 $\lambda_i > 0$ | 特征值全正 |
| 正定 $\Leftrightarrow$ Sylvester 准则 | 顺序主子式全正 |
| 正定 $\Leftrightarrow$ Cholesky 存在 | $A = LL^T$，高效数值方法 |
| 能量函数 | $E = \frac{1}{2}\mathbf{x}^TK\mathbf{x} > 0$ → 稳定 |
| Lyapunov 稳定性 | $A^TP + PA < 0$，$P > 0$ → 系统稳定 |
| 凸优化 | Hessian 正定 → 严格凸 → 唯一极小值 |
| 条件数 $\kappa$ | $\kappa$ 越大，优化越难，收敛越慢 |
| 封闭性 | 正定矩阵的和、逆、幂等仍正定 |

---

## 习题

### 概念理解

**7.1** 判断正误：
  - (a) 所有对角元素为正的对称矩阵一定正定。
  - (b) 正定矩阵的行列式一定为正。
  - (c) 若 $A$ 正定，则 $A + I$ 也正定。
  - (d) 两个正定矩阵的乘积 $AB$ 一定正定。

**7.2** 解释：为什么 Cholesky 分解只适用于正定矩阵？如果矩阵半正定会发生什么？

### 计算练习

**7.3** 对以下矩阵判别正定性（用两种方法），若正定则求 Cholesky 分解：

$$
(a) \quad A = \begin{pmatrix} 9 & 6 \\ 6 & 5 \end{pmatrix} \qquad
(b) \quad B = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 3 & 1 \\ 0 & 1 & 2 \end{pmatrix}
$$

**7.4** 求矩阵 $A = \begin{pmatrix} 4 & 2 \\ 2 & 3 \end{pmatrix}$ 的正定平方根 $A^{1/2}$。验证 $A^{1/2}A^{1/2} = A$。

**7.5** 设 $A = \begin{pmatrix} -3 & 1 \\ -1 & -2 \end{pmatrix}$。利用 Lyapunov 方程 $A^TP + PA = -I$ 求 $P$，并验证 $P$ 正定。

### 思考题

**7.6** 正定矩阵 $A$ 和 $B$ 的乘积 $AB$ 不一定正定（因为 $AB$ 不一定对称）。但如果 $AB = BA$，证明 $AB$ 正定。

**7.7** 在机器学习中，为什么正则化项 $\alpha\|\mathbf{w}\|^2$ 能改善优化问题的条件数？用 Hessian 矩阵说明。

### 编程题

**7.8** 数值比较三种方法求解正定线性方程组 $A\mathbf{x} = \mathbf{b}$ 的效率：
  - 直接 LU 分解 (`numpy.linalg.solve`)
  - Cholesky 分解 (`scipy.linalg.cho_solve`)
  - 迭代法（共轭梯度 `scipy.sparse.linalg.cg`）
  - 对 $n = 100, 500, 1000$ 的正定矩阵进行测试，记录时间和误差

**7.9** 实现条件数对梯度下降影响的完整实验：
  - 构造条件数分别为 1, 10, 100, 1000 的正定矩阵
  - 对每种条件数，运行梯度下降直到收敛
  - 绘制：(a) 收敛步数 vs 条件数，(b) 各种条件数下的等高线与迭代轨迹

---

> **下一章预告**：分块矩阵将大问题拆分为小问题——每个小问题可以独立处理。块对角化是"分而治之"思想在线性代数中的体现。
