
## A. $2 \times 2$ 矩阵对角化

考虑一个二能级系统，其哈密顿量为：

$$H = \begin{bmatrix} a & c \\ c^* & -a \end{bmatrix},$$

其中 $a$ 是实数，$c$ 是复数。

1. 对角化哈密顿量，作分解 $H = Q \Lambda Q^{-1}$。其中 $\Lambda$ 为特征值构成的对角矩阵，$Q$ 为本征列向量构成的矩阵，且所有本征向量归一化。$Q^{-1}$ 和 $Q$ 有什么关系，为什么？（1分）
    
2. 证明公式：$e^{iHt} = Q e^{i\Lambda t} Q^{-1}$，并用这一公式计算出 $e^{iHt}$ 的具体表达式。（1分）
    

---

## B. 横场 Ising 模型

研究一个自旋1/2的横场Ising模型：

$$H = -J \sum_{i=1}^{L} \hat{\sigma}_i^z \hat{\sigma}_{i+1}^z - h \sum_{i=1}^{L} \hat{\sigma}_i^x.$$

取周期边界条件（$\hat{\sigma}_{L+1}^z \equiv \hat{\sigma}_1^z$）且固定 $J = 1$。

1. 取参数 $L = 12$，对于不同的横场强度 $h = 0.5, 1.0, 2.0$，分别对角化哈密顿量 $H$，求得**基态**和**第一激发态**并列出其能量（2分）。对于基态波函数，计算第一个格点上横场方向磁化强度的期望 $\langle \sigma_1^x \rangle$，并一并列出。（1分）
    
2. $t = 0$ 时系统处于 $h = 0.5$ 的基态，$t = 0^+$ 的瞬间哈密顿量参数变为 $h = 3.0$，求在此哈密顿量下的时间演化，计算 $\langle \sigma_1^x(t) \rangle$ 在时间 $t \in [0, 20]$ 的变化情况，并画出来。你需要给出精确的结果和 Runge-Kutta 方法的结果。验证两种方法的结果是一样的。（3分）
    
3. 取参数 $L = 18, h = 1$。从Néel态 $|\psi(0)\rangle = |\uparrow \downarrow \uparrow \downarrow \cdots\rangle$（奇数格点自旋朝上，偶数格点朝下）出发，求解波函数 $|\psi(t)\rangle$ 在区间 $t \in [0, 10]$ 上面的时间演化。画出 Fidelity $F(t) = |\langle \psi(0) | \psi(t) \rangle|^2$ 随时间的变化。（2分）