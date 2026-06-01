# 计算物理导论 - Homework 10 : 全连接Ising模型

考虑一个 $N$ 个格点的全连接Ising模型，其哈密顿量为：

$$
H = -\frac{1}{2}\sum_{i\neq j} J_{ij}s_i s_j,
$$

其中 $s_i = \pm 1$ 为格点的自旋状态，耦合强度矩阵 $J_{ij}$ 是对称的。接下来让我们定义 $p$ 个自旋态集合：

$$
\{ \mathbf{v}^\mu \} = \{ \mathbf{v}^1, \mathbf{v}^2, \cdots, \mathbf{v}^p \},
$$

$\mathbf{v}^\mu$ 是可以任意选取的 $N$ 格点的自旋组态。对于一个给定的集合 $\{ \mathbf{v}^\mu \}$，按照下述规则确定一个耦合强度矩阵：

$$
\left. J_{ij} \right|_{i\neq j} = \frac{1}{N}\sum_{\mu=1}^p v_i^\mu v_j^\mu, \quad J_{ii} = 0.
$$

定义两个自旋构型 $\mathbf{u, v}$ 的余弦相似度：

$$
S_C(\mathbf{u, v}) = \frac{\left|\sum_i^N u_i v_i\right|}{\sqrt{\sum_i^N u_i^2} \cdot \sqrt{\sum_i^N v_i^2}} = \frac{1}{N} \left|\sum_i^N u_i v_i\right|.
$$

1. 考虑最简单的情况 $p = 1$，即初始的集合 $\{ \mathbf{v}^\mu \}$ 只有一个元素 $\mathbf{v} = [v_1, v_2, \cdots, v_N]^T$。写出一个任意自旋态 $\mathbf{s} = [s_1, s_2, \cdots, s_N]^T$ 对应的能量 $E(\mathbf{s})$。系统的基态是什么？基态简并吗？（2分）
2. 产生一个大小为 $p$ 的集合 $\{ \mathbf{v}^\mu \}$。其中 $v_i^\mu = \text{rand}(-1, +1)$ 随机选取。
    i. 取值 $N = 1000, p = 1$，且随机生成一个初态 $\mathbf{s}$，它和集合里唯一元素的相似度 $S_C(\mathbf{s, v})$ 是多少？现在用metropolis单自旋更新的方法来模拟 $H$ 的基态（$\beta \rightarrow \infty$）。当系统能量不再变化的时候，你找到的 $\mathbf{s}$ 是什么，有什么特点？ $S_C(\mathbf{s, v})$ 又是多少？（2分）
    ii. 取值 $N = 2000, p = 10$。产生一个初态集合 $\{ \mathbf{s}^\mu \}$，其中 $\mathbf{s}^\mu$ 为将 $\mathbf{v}^\mu$ 加上噪声得到，其伪代码为：
    
    ```text
    add_noise(v, δ) = [(rand()<δ ? rand((-1,+1)) : vi) for vi in v]
    ```

    其中 $\delta \in [0, 1]$ 是噪声强度。现在取 $\delta = 0.5$。对得到的集合 $\{ \mathbf{s}^\mu \}$ 都作演化，给出末态和对应初态的相似度 $S_C^\mu = S_C(\mathbf{u}^\mu, \mathbf{v}^\mu)$。你发现了什么？（2分）
    
    iii. 保持 $N = 2000$ 并改变 $p$ 的值。对于不同的 $p$ 值，计算 $S_C^\mu$ 的分布，并用直方图画出来。这一分布随着 $p$ 的值改变出现了什么变化？这样的变化说明了什么？（2分） hint：$p$ 在此题可以在 $50 - 400$ 间选取。为了得到更好的分布，你可以考虑对不同的 $\{ \mathbf{v}^\mu \}$ 作系统平均。
    
    iv. 你能解释你的发现吗？（2分） 你还有什么别的思考？