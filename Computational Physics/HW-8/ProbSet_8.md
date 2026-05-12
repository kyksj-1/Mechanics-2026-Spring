# 计算物理导论 - Homework 8: Monte-Carlo

**注：** 本作业中默认取无量纲单位 $k_B = 1, T = 1/\beta$。

## A. Ising 模型

使用 Monte Carlo 方法模拟 $L \times L$ 二维方晶格上的经典 Ising 模型：

$$H = -\sum_{\langle ij \rangle} J_{ij} \sigma_i \sigma_j$$

其中 $\langle ij \rangle$ 取不重复的最近邻邻居，且固定 $J_{ij} = J = 1$。对晶格取周期性边界条件。

1. **精确计算：** 对 $L = 4$ 的情况，精确计算 $T = 1$ 情况下的平衡态能量 $E$ 和自由能 $F$。（1分）
    
2. **理论推导：** 写出一般情况 Markov Chain Monte Carlo (MCMC) 的细致平衡方程。对于 Ising 模型，构型的权重是什么？你选择的更新方法有哪些过程（过程和逆过程）？根据过程两侧的状态权重，设计一个选择概率，并计算接受概率。（1分）
    
3. **数值模拟：** 使用 Monte-Carlo 计算 $L = 4, T = 1$ 的平衡态能量 $\langle E \rangle$。验证你的算法是正确的。（1分）
    
4. **相变观测：** 计算 $L = 8, 16, 32$ 的物理量随温度变化的关系。温度区间取 $[1.5, 3]$，间距为 $0.1$。要计算的物理量包括：
    
    - 磁化强度平方 $\langle m^2 \rangle = \frac{\langle M^2 \rangle}{N^2}$
        
    - 比热 $c = \frac{\beta^2 (\langle E^2 \rangle - \langle E \rangle^2)}{N}$
        
    - 磁化率 $\chi = \frac{\beta (\langle M^2 \rangle - \langle |M| \rangle^2)}{N}$
        
    
    对每个物理量，将不同 $L$ 的结果画在同一张图中。你发现了什么？（2分）
    

---

## B. 弛豫动力学

仍然考虑 (A) 中的模型，固定更新算法为：

- 每次更新在晶格中随机选取一个格点，尝试进行标准的 Metropolis 更新。
    
- 每随机尝试更新 $L^2$ 次定义为一个蒙卡步（MCS）。
    

初始化无穷高温的系统，并骤临界温度：

$$T_c = \frac{2}{\ln(1 + \sqrt{2})}$$

进行演化。计算系统的平均能量 $\langle E(t) \rangle$，其中 $t$ 是蒙卡时间步。

1. **弛豫过程：** 对 $L = 16$ 的系统，画出能量随时间的变化关系。粗略探究需要多长时间，系统能量弛豫到稳态 $\langle E(\infty) \rangle$。（2分）
    
2. **尺度效应：** 改变系统的尺寸，观察系统能量相对稳态的差距 $\Delta(t) = \langle E(t) \rangle - \langle E(\infty) \rangle$ 的长时间行为。你发现了什么规律？系统尺寸对这个规律有怎样的影响？临界温度在这个问题中可能有什么意义？（3分）
    
    - **Hint:** 谨慎地确定 $\langle E(\infty) \rangle$。