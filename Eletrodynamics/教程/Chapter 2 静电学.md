

# 第2章 静电学 (Electrostatics)

## 引言：从力到场的范式转移

在经典力学中，我们习惯于"超距作用"的思维：两个质量之间的引力瞬间传递，不需要任何中间媒介。库仑在1785年发现电荷之间的力也服从类似的平方反比律时，人们自然地沿用了同样的思路。

但这幅图景是有问题的。如果一个电荷突然移动，远处的电荷是否**瞬间**感受到力的变化？答案是否定的——信息的传播速度不能超过光速。那么，力是如何从一个电荷"传递"到另一个电荷的？

答案是**场（Field）**。电荷 $q$ 不直接对远处的电荷施力，而是在自己周围的空间中创造了一种物理实在——**电场** $\mathbf{E}$。其他电荷感受到的力，来自它所在位置的电场，而非来自遥远的源电荷本身。

本章处理的是**静电学（Electrostatics）**：所有电荷都静止不动的情形。在这个简化但基本的场景中，我们将：
1. 建立电场的定义和计算方法；
2. 发现电场的两个微分性质（散度与旋度），得到高斯定律；
3. 引入电势，将矢量问题化为标量问题；
4. 探讨电场储存的能量；
5. 理解导体在电场中的行为。

这些结果构成了整个电动力学大厦的地基。

---

## 2.1 电场

### 2.1.1 库仑定律：一切的起点

实验事实表明，两个静止点电荷 $q_1$ 和 $q_2$ 之间的力为：

$$\mathbf{F} = \frac{1}{4\pi\varepsilon_0} \frac{q_1 q_2}{\mathcal{r}^2} \hat{\boldsymbol{\mathcal{r}}}$$

其中 $\boldsymbol{\mathcal{r}} = \mathbf{r}_1 - \mathbf{r}_2$ 是从 $q_2$ 指向 $q_1$ 的分离矢量（回忆第1章 §1.1.4 的记号），$\varepsilon_0 = 8.854 \times 10^{-12} \text{ C}^2/(\text{N}\cdot\text{m}^2)$ 是真空电容率。

**库仑定律的三个关键特征**：

1. **平方反比**：力与距离的平方成反比。这不是巧合——它与三维空间的几何结构深刻相关（后续高斯定律将揭示这一点）。
2. **线性叠加**：多个电荷对某一电荷的力是各自力的矢量和，**不受其他电荷的存在影响**。这是一个非凡的实验事实，绝非理所当然。
3. **与 $1/r^2$ 的精确偏差**：实验已将库仑定律的指数验证到 $|2 - n| < 10^{-16}$。任何偏差都将意味着光子具有非零质量。

### 2.1.2 电场的定义与叠加原理

我们引入**电场**的概念。设源电荷 $q$ 位于 $\mathbf{r}'$，则它在空间中场点 $\mathbf{r}$ 处产生的电场为：

$$\mathbf{E}(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0} \frac{q}{\mathcal{r}^2} \hat{\boldsymbol{\mathcal{r}}}$$

其中 $\boldsymbol{\mathcal{r}} = \mathbf{r} - \mathbf{r}'$。

**物理解读**：电场 $\mathbf{E}(\mathbf{r})$ 是源电荷在空间每一点留下的"指令"——如果一个试探电荷 $Q$ 被放到 $\mathbf{r}$ 处，它将感受到力 $\mathbf{F} = Q\mathbf{E}(\mathbf{r})$。

对于 $N$ 个点电荷 $q_1, q_2, \dots, q_N$ 分别位于 $\mathbf{r}_1', \mathbf{r}_2', \dots, \mathbf{r}_N'$，**叠加原理**给出总电场：

$$\mathbf{E}(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0} \sum_{i=1}^{N} \frac{q_i}{\mathcal{r}_i^2} \hat{\boldsymbol{\mathcal{r}}}_i$$

其中 $\boldsymbol{\mathcal{r}}_i = \mathbf{r} - \mathbf{r}_i'$。

**关于叠加原理的深层意义**：叠加原理意味着麦克斯韦方程组是**线性**的。这使得我们可以将任何复杂的电荷分布分解为简单部分分别求解，然后叠加。非线性介质（如铁电体）中叠加原理不再成立，问题会变得极其困难。

### 2.1.3 连续电荷分布

当电荷数量极大（如一根带电导线上有 $\sim 10^{23}$ 个多余电子）时，离散求和变为积分。根据电荷的几何分布，我们有三种情况：

**线电荷密度 $\lambda$（单位：C/m）**：电荷沿一维曲线分布。

$$\mathbf{E}(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0} \int \frac{\lambda(\mathbf{r}')}{\mathcal{r}^2} \hat{\boldsymbol{\mathcal{r}}} \, dl'$$

**面电荷密度 $\sigma$（单位：C/m²）**：电荷在二维曲面上分布。

$$\mathbf{E}(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0} \int \frac{\sigma(\mathbf{r}')}{\mathcal{r}^2} \hat{\boldsymbol{\mathcal{r}}} \, da'$$

**体电荷密度 $\rho$（单位：C/m³）**：电荷在三维区域内分布。

$$\mathbf{E}(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0} \int \frac{\rho(\mathbf{r}')}{\mathcal{r}^2} \hat{\boldsymbol{\mathcal{r}}} \, d\tau'$$

**计算警示**：这些积分中的 $\hat{\boldsymbol{\mathcal{r}}}$ 和 $\mathcal{r}$ 都是积分变量 $\mathbf{r}'$ 的函数，不能提到积分号外。此外，$\hat{\boldsymbol{\mathcal{r}}}$ 是矢量，意味着我们通常需要分量来做积分。

---

**例题 2.1**：求距均匀带电直线段（长度 $2L$，线电荷密度 $\lambda$）中点正上方距离 $z$ 处的电场。

**解**：建立坐标系，令线段沿 $x$ 轴放置，中点在原点，场点 $P$ 在 $z$ 轴上 $(0, 0, z)$。

在位置 $x$ 处取微元 $dx$，该微元到 $P$ 的距离 $\mathcal{r} = \sqrt{x^2 + z^2}$。由对称性，$x$ 处和 $-x$ 处微元产生的电场水平分量相消，只剩 $z$ 分量：

$$dE_z = \frac{1}{4\pi\varepsilon_0} \frac{\lambda \, dx}{x^2 + z^2} \cdot \frac{z}{\sqrt{x^2 + z^2}}$$

积分：

$$E_z = \frac{\lambda z}{4\pi\varepsilon_0} \int_{-L}^{L} \frac{dx}{(x^2 + z^2)^{3/2}} = \frac{\lambda z}{4\pi\varepsilon_0} \cdot \frac{2L}{z^2\sqrt{L^2 + z^2}}$$

$$\boxed{E = \frac{1}{4\pi\varepsilon_0} \frac{2\lambda L}{z\sqrt{z^2 + L^2}} \hat{\mathbf{z}}}$$

**验证极限**：
- $z \gg L$（远处）：$E \approx \frac{1}{4\pi\varepsilon_0}\frac{2\lambda L}{z^2} = \frac{1}{4\pi\varepsilon_0}\frac{q}{z^2}$，其中 $q = 2\lambda L$ 是总电荷。远处线段看起来像点电荷。
- $L \to \infty$（无限长线）：$E = \frac{\lambda}{2\pi\varepsilon_0 z}$，只与距离的一次方成反比——这是无限长直线的经典结果。

---

**例题 2.2**：求均匀带电圆环（半径 $R$，总电荷 $Q$）轴线上距中心 $z$ 处的电场。

**解**：设圆环位于 $xy$ 平面，圆心在原点。由对称性，轴线上只有 $z$ 分量存在。

圆环上微元 $dl' = R \, d\phi'$，其到场点的距离为 $\mathcal{r} = \sqrt{R^2 + z^2}$（对所有微元相同）。

$$E_z = \frac{1}{4\pi\varepsilon_0} \int_0^{2\pi} \frac{\lambda R \, d\phi'}{R^2 + z^2} \cdot \frac{z}{\sqrt{R^2 + z^2}} = \frac{1}{4\pi\varepsilon_0} \frac{Qz}{(R^2 + z^2)^{3/2}}$$

$$\boxed{\mathbf{E} = \frac{1}{4\pi\varepsilon_0} \frac{Qz}{(R^2 + z^2)^{3/2}} \hat{\mathbf{z}}}$$

**注意**：$z = 0$ 时 $E = 0$（圆心处电场为零，因为对称性使所有贡献相消）；$z \gg R$ 时 $E \approx Q/(4\pi\varepsilon_0 z^2)$（回到点电荷）。

---

**例题 2.3**：求均匀带电圆盘（半径 $R$，面电荷密度 $\sigma$）轴线上距中心 $z$ 处的电场。

**解**：将圆盘视为无穷多同心圆环的叠加。半径 $r'$、宽度 $dr'$ 的细圆环携带电荷 $dq = \sigma \cdot 2\pi r' dr'$。利用例题 2.2 的结果：

$$E_z = \frac{\sigma z}{4\pi\varepsilon_0} \int_0^R \frac{2\pi r' \, dr'}{(r'^2 + z^2)^{3/2}}$$

令 $u = r'^2 + z^2$，$du = 2r' dr'$：

$$E_z = \frac{\sigma z}{4\pi\varepsilon_0} \cdot \pi \int_{z^2}^{R^2+z^2} u^{-3/2} du = \frac{\sigma z}{4\pi\varepsilon_0} \cdot 2\pi \left[\frac{1}{|z|} - \frac{1}{\sqrt{R^2+z^2}}\right]$$

对于 $z > 0$：

$$\boxed{E_z = \frac{\sigma}{2\varepsilon_0}\left[1 - \frac{z}{\sqrt{R^2 + z^2}}\right]}$$

**极限验证**：
- $R \to \infty$（无限大平面）：$E_z = \sigma/(2\varepsilon_0)$，与距离无关！这是无限大平面的经典结果。
- $z \gg R$：展开 $1/\sqrt{1 + R^2/z^2} \approx 1 - R^2/(2z^2)$，得 $E \approx \sigma R^2/(4\varepsilon_0 z^2) = Q/(4\pi\varepsilon_0 z^2)$，回到点电荷。

---

## 2.2 静电场的散度与旋度

原则上，§2.1 的积分公式解决了静电学的全部问题——知道电荷分布就能算出电场。但这些积分通常很难计算。幸运的是，电场的**微分性质**（散度和旋度）提供了更优雅的工具。

### 2.2.1 电场线与电通量

**电场线**是这样的曲线族：在每一点上，曲线的切线方向与电场方向一致。它们从正电荷出发，终止于负电荷（或延伸到无穷远）。

电场线的**密度**（单位面积穿过的条数）正比于场强 $|\mathbf{E}|$。因此，对于点电荷 $q$，从 $q$ 发出的场线均匀辐射，穿过半径为 $r$ 的球面的条数密度正比于 $q/r^2$，恰好与 $|\mathbf{E}| \propto 1/r^2$ 一致。

**电通量**（Electric Flux）定义为电场通过某个曲面的积分：

$$\Phi_E = \int_S \mathbf{E} \cdot d\mathbf{a}$$

对于闭合曲面，通量的物理意义是：穿出曲面的场线净条数。

### 2.2.2 高斯定律的推导

**直觉**：考虑一个以点电荷 $q$ 为中心、半径为 $r$ 的球面。穿过该球面的电通量为：

$$\oint_S \mathbf{E} \cdot d\mathbf{a} = \int \frac{1}{4\pi\varepsilon_0}\frac{q}{r^2}\hat{\mathbf{r}} \cdot r^2\sin\theta \, d\theta \, d\phi \, \hat{\mathbf{r}} = \frac{q}{4\pi\varepsilon_0} \cdot 4\pi = \frac{q}{\varepsilon_0}$$

$r^2$ 的消去是关键——面积增长 $\propto r^2$ 恰好补偿了场强衰减 $\propto 1/r^2$。这也是为什么高斯定律只在三维空间的平方反比力中成立。

**对任意闭合曲面的推广**：

关键论证如下：
1. 对于以 $q$ 为中心的球面，$\oint \mathbf{E} \cdot d\mathbf{a} = q/\varepsilon_0$。
2. 对于包含 $q$ 的任意闭合曲面，每条电场线既然从 $q$ 出发，必须穿出曲面。无论曲面形状如何，穿出的净条数不变。（严格证明利用立体角：$\mathbf{E} \cdot d\mathbf{a} = (q/4\pi\varepsilon_0) d\Omega$，而闭合曲面对内部点张的总立体角为 $4\pi$。）
3. 若 $q$ 在曲面外部，从 $q$ 发出的场线穿入曲面后必定穿出，净贡献为零。
4. 由叠加原理，对多个电荷：

$$\oint_S \mathbf{E} \cdot d\mathbf{a} = \frac{1}{\varepsilon_0} \sum_{i \in \text{内部}} q_i = \frac{Q_{\text{enc}}}{\varepsilon_0}$$

这就是**高斯定律（积分形式）**：

$$\boxed{\oint_S \mathbf{E} \cdot d\mathbf{a} = \frac{Q_{\text{enc}}}{\varepsilon_0}}$$

**高斯定律（微分形式）**：

将 $Q_{\text{enc}} = \int_V \rho \, d\tau$ 代入，并对左边使用散度定理：

$$\oint_S \mathbf{E} \cdot d\mathbf{a} = \int_V (\nabla \cdot \mathbf{E}) \, d\tau = \frac{1}{\varepsilon_0}\int_V \rho \, d\tau$$

由于这对**任意**体积 $V$ 成立，被积函数必须处处相等：

$$\boxed{\nabla \cdot \mathbf{E} = \frac{\rho}{\varepsilon_0}}$$

这是麦克斯韦方程组的第一个方程。它的物理含义极其清晰：**电荷是电场的源。正电荷是源（场线涌出），负电荷是汇（场线汇入）。**

**数学验证**：对点电荷 $\mathbf{E} = \frac{q}{4\pi\varepsilon_0}\frac{\hat{\mathbf{r}}}{r^2}$，利用第1章的结论 $\nabla \cdot (\hat{\mathbf{r}}/r^2) = 4\pi\delta^3(\mathbf{r})$，得到 $\nabla \cdot \mathbf{E} = (q/\varepsilon_0)\delta^3(\mathbf{r})$，与 $\rho = q\delta^3(\mathbf{r})$ 完全一致。

### 2.2.3 高斯定律的应用

高斯定律 $\oint \mathbf{E} \cdot d\mathbf{a} = Q_{\text{enc}}/\varepsilon_0$ 本身是精确的，但作为**计算工具**，它只在具有高度对称性的问题中有效——此时我们可以选择适当的高斯面，使 $\mathbf{E}$ 在面上为常数或为零，从而将积分化为代数。

三种经典对称性：

---

**（1）球对称**：$\rho$ 只依赖于 $r$

选择以对称中心为球心、半径为 $r$ 的球面作为高斯面。由对称性，$\mathbf{E} = E(r)\hat{\mathbf{r}}$，在球面上 $E$ 为常数。

$$\oint \mathbf{E} \cdot d\mathbf{a} = E(r) \cdot 4\pi r^2 = \frac{Q_{\text{enc}}(r)}{\varepsilon_0}$$

$$\boxed{E(r) = \frac{Q_{\text{enc}}(r)}{4\pi\varepsilon_0 r^2}}$$

**例题 2.4**：均匀带电实心球（半径 $R$，体电荷密度 $\rho_0$，总电荷 $Q = \frac{4}{3}\pi R^3 \rho_0$）。

- 球外（$r > R$）：$Q_{\text{enc}} = Q$，故 $E = \frac{Q}{4\pi\varepsilon_0 r^2}$，与点电荷无异。
- 球内（$r < R$）：$Q_{\text{enc}} = \rho_0 \cdot \frac{4}{3}\pi r^3 = Q\frac{r^3}{R^3}$，故 $E = \frac{Qr}{4\pi\varepsilon_0 R^3}$，线性增长。

$$\boxed{E(r) = \begin{cases} \dfrac{Qr}{4\pi\varepsilon_0 R^3} & r \leq R \\[6pt] \dfrac{Q}{4\pi\varepsilon_0 r^2} & r \geq R \end{cases}}$$

**物理图景**：球内只有半径 $r$ 以内的电荷对该处电场有贡献，外层电荷的贡献**精确为零**（球壳定理）。

---

**（2）柱对称**：$\rho$ 只依赖于 $s$（到轴线距离）

选择同轴圆柱面（半径 $s$，高度 $h$）作为高斯面。由对称性，$\mathbf{E} = E(s)\hat{\mathbf{s}}$。

$$\oint \mathbf{E} \cdot d\mathbf{a} = E(s) \cdot 2\pi s h = \frac{Q_{\text{enc}}}{\varepsilon_0}$$

**例题 2.5**：无限长均匀带电直线（线电荷密度 $\lambda$）。

$$E(s) = \frac{\lambda}{2\pi\varepsilon_0 s}$$

$$\boxed{\mathbf{E} = \frac{\lambda}{2\pi\varepsilon_0 s}\hat{\mathbf{s}}}$$

验证：这与例题 2.1 中 $L \to \infty$ 的极限一致。但高斯定律在几行之内就得到了答案，而直接积分需要繁琐得多的计算。

---

**（3）平面对称**：$\rho$ 只依赖于某一坐标（如 $z$）

选择"高斯药盒"（Gaussian pillbox）——一个扁平圆柱，上下底面平行于带电平面。

**例题 2.6**：无限大均匀带电平面（面电荷密度 $\sigma$）。

由对称性，$\mathbf{E}$ 垂直于平面，上方向上，下方向下。药盒的侧面对通量无贡献，两个底面（面积 $A$）各贡献 $EA$：

$$2EA = \frac{\sigma A}{\varepsilon_0}$$

$$\boxed{\mathbf{E} = \frac{\sigma}{2\varepsilon_0}\hat{\mathbf{n}}}$$

其中 $\hat{\mathbf{n}}$ 是指向远离平面一侧的法向量。

**惊人之处**：电场与到平面的距离无关！这是无限大平面的特性。在现实中，只有在距平面的距离远小于平面尺寸时，此结果才是好近似。

---

**例题 2.7**：两块平行的无限大平面，分别携带 $+\sigma$ 和 $-\sigma$。

由叠加，两板之间场强加倍，两板外侧相消：

$$\mathbf{E} = \begin{cases} 0 & \text{两板外侧} \\ \sigma/\varepsilon_0 & \text{两板之间} \end{cases}$$

这就是**平行板电容器**的理想化模型。

### 2.2.4 静电场的旋度

我们来计算静电场的旋度。从点电荷的电场出发：

$$\mathbf{E} = \frac{q}{4\pi\varepsilon_0}\frac{\hat{\mathbf{r}}}{r^2}$$

回忆第1章习题 [G 1.13] 的结论：$\hat{\mathbf{r}}/r^2 = -\nabla(1/r)$。因此 $\mathbf{E} = -\frac{q}{4\pi\varepsilon_0}\nabla(1/r)$，这是一个梯度场。

由矢量恒等式 $\nabla \times (\nabla f) = 0$（任何梯度场无旋），得到：

$$\nabla \times \mathbf{E} = 0$$

由叠加原理，多个电荷的总电场的旋度也为零。对连续分布同理。

**积分形式**（由斯托克斯定理）：

$$\boxed{\oint_C \mathbf{E} \cdot d\mathbf{l} = 0}$$

即静电场沿任意闭合回路的环流为零。这意味着静电力是**保守力**，做功只与起点和终点有关，与路径无关。

**物理意义**：如果沿某闭合路径环流不为零，就可以让电荷不断绕圈获取能量——永动机！静电场的无旋性正是能量守恒的体现。

**重要注记**：$\nabla \times \mathbf{E} = 0$ **只对静电场成立**。在时变情况下（第7章），变化的磁场会产生有旋的感生电场（法拉第定律），此时 $\nabla \times \mathbf{E} = -\partial\mathbf{B}/\partial t \neq 0$。

---

## 2.3 电势

### 2.3.1 从无旋性到标量势

既然 $\nabla \times \mathbf{E} = 0$，由亥姆霍兹定理（第1章 §1.6），$\mathbf{E}$ 可以写为某个标量场的梯度：

$$\boxed{\mathbf{E} = -\nabla V}$$

标量函数 $V(\mathbf{r})$ 称为**电势（Electric Potential）**。负号是约定：电场从高电势指向低电势（正电荷从"山顶"滚向"山谷"）。

**为什么电势如此有用？**

1. **从三个分量到一个函数**：$\mathbf{E}$ 有三个分量 $(E_x, E_y, E_z)$，而 $V$ 只是一个标量。知道 $V$ 后，取梯度即可恢复 $\mathbf{E}$ 的全部信息。
2. **积分更简单**：计算电势时，被积函数是标量 $1/\mathcal{r}$，不需要处理烦人的矢量分量。
3. **能量的直接体现**：电势差就是单位电荷的势能差。

### 2.3.2 电势与电场的关系

**从电场到电势**（积分关系）：

由梯度定理 $\int_a^b (\nabla V) \cdot d\mathbf{l} = V(b) - V(a)$，代入 $\mathbf{E} = -\nabla V$：

$$V(\mathbf{r}) - V(\mathbf{r}_{\text{ref}}) = -\int_{\mathbf{r}_{\text{ref}}}^{\mathbf{r}} \mathbf{E} \cdot d\mathbf{l}$$

选取参考点 $V(\mathbf{r}_{\text{ref}}) = 0$（对有限电荷分布，通常取无穷远为零势点）：

$$\boxed{V(\mathbf{r}) = -\int_\infty^{\mathbf{r}} \mathbf{E} \cdot d\mathbf{l}}$$

**注意**：由于 $\nabla \times \mathbf{E} = 0$，积分与路径无关，$V$ 是良定义的。

**两点之间的电势差**：

$$V(\mathbf{b}) - V(\mathbf{a}) = -\int_{\mathbf{a}}^{\mathbf{b}} \mathbf{E} \cdot d\mathbf{l}$$

这是电压表测量的物理量——单位正电荷从 $\mathbf{a}$ 移到 $\mathbf{b}$ 时，电场力做的功为 $q[V(\mathbf{a}) - V(\mathbf{b})]$。

### 2.3.3 电势的计算

**（1）由电场积分求电势**

若已知 $\mathbf{E}$，直接使用 $V(\mathbf{r}) = -\int_\infty^{\mathbf{r}} \mathbf{E} \cdot d\mathbf{l}$。

**例题 2.8**：均匀带电球壳（半径 $R$，总电荷 $Q$）的电势。

球壳外 $\mathbf{E} = \frac{Q}{4\pi\varepsilon_0 r^2}\hat{\mathbf{r}}$，球壳内 $\mathbf{E} = 0$。

球壳外（$r > R$）：
$$V(r) = -\int_\infty^r \frac{Q}{4\pi\varepsilon_0 r'^2} dr' = \frac{Q}{4\pi\varepsilon_0 r}$$

球壳内（$r < R$）：从无穷远积分到 $R$，再从 $R$ 积分到 $r$（后一段 $\mathbf{E} = 0$）：
$$V(r) = -\int_\infty^R \frac{Q}{4\pi\varepsilon_0 r'^2}dr' - \int_R^r 0 \, dr' = \frac{Q}{4\pi\varepsilon_0 R}$$

$$\boxed{V(r) = \begin{cases} \dfrac{Q}{4\pi\varepsilon_0 R} & r \leq R \\[6pt] \dfrac{Q}{4\pi\varepsilon_0 r} & r \geq R \end{cases}}$$

球壳内电势为常数（等于球面上的电势），尽管内部电场为零。

---

**（2）由电荷分布直接计算电势**

对点电荷 $q$ 位于原点：$V(r) = \frac{q}{4\pi\varepsilon_0 r}$。

对一般电荷分布，利用叠加：

$$\boxed{V(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0} \int \frac{\rho(\mathbf{r}')}{\mathcal{r}} d\tau'}$$

类似地有 $V = \frac{1}{4\pi\varepsilon_0}\int\frac{\lambda}{\mathcal{r}}dl'$ 和 $V = \frac{1}{4\pi\varepsilon_0}\int\frac{\sigma}{\mathcal{r}}da'$。

**与电场公式的对比**：电势公式中是 $1/\mathcal{r}$（标量），而电场公式中是 $\hat{\boldsymbol{\mathcal{r}}}/\mathcal{r}^2$（矢量）。这就是电势计算更简便的原因。

---

**例题 2.9**：均匀带电实心球（半径 $R$，总电荷 $Q$）的电势。

球外：$V(r) = \frac{Q}{4\pi\varepsilon_0 r}$（等效于点电荷）。

球内：从无穷远积分到 $r$：

$$V(r) = -\int_\infty^R \frac{Q}{4\pi\varepsilon_0 r'^2}dr' - \int_R^r \frac{Q r'}{4\pi\varepsilon_0 R^3}dr'$$

$$= \frac{Q}{4\pi\varepsilon_0 R} - \frac{Q}{4\pi\varepsilon_0 R^3}\cdot\frac{r^2 - R^2}{2} = \frac{Q}{4\pi\varepsilon_0}\cdot\frac{3R^2 - r^2}{2R^3}$$

$$\boxed{V(r) = \begin{cases} \dfrac{Q}{4\pi\varepsilon_0}\dfrac{3R^2 - r^2}{2R^3} & r \leq R \\[6pt] \dfrac{Q}{4\pi\varepsilon_0 r} & r \geq R \end{cases}}$$

**验证**：$\mathbf{E} = -\nabla V = -\frac{dV}{dr}\hat{\mathbf{r}}$。球内：$-\frac{d}{dr}\left[\frac{Q(3R^2 - r^2)}{8\pi\varepsilon_0 R^3}\right] = \frac{Qr}{4\pi\varepsilon_0 R^3}$，与例题 2.4 一致。

### 2.3.4 泊松方程与拉普拉斯方程

将 $\mathbf{E} = -\nabla V$ 代入高斯定律 $\nabla \cdot \mathbf{E} = \rho/\varepsilon_0$：

$$\nabla \cdot (-\nabla V) = \frac{\rho}{\varepsilon_0}$$

$$\boxed{\nabla^2 V = -\frac{\rho}{\varepsilon_0}} \qquad \text{（泊松方程）}$$

在无电荷区域（$\rho = 0$），退化为：

$$\boxed{\nabla^2 V = 0} \qquad \text{（拉普拉斯方程）}$$

泊松方程是电势理论的核心方程。它的解是：

$$V(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0}\int \frac{\rho(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} d\tau'$$

这可以通过 $\nabla^2(1/|\mathbf{r} - \mathbf{r}'|) = -4\pi\delta^3(\mathbf{r} - \mathbf{r}')$（第1章 §1.5.4 的推广）直接验证。

### 2.3.5 边界条件

当电场穿越带有面电荷密度 $\sigma$ 的界面时，电场不连续。利用高斯定律和环路定理可以精确确定跃变量。

**法向分量的跃变**：

在界面两侧取高斯药盒（厚度 $\to 0$），仅上下底面贡献：

$$E_{\text{above}}^{\perp} \cdot A - E_{\text{below}}^{\perp} \cdot A = \frac{\sigma A}{\varepsilon_0}$$

$$\boxed{E_{\text{above}}^{\perp} - E_{\text{below}}^{\perp} = \frac{\sigma}{\varepsilon_0}}$$

**切向分量的连续性**：

沿界面取一个窄矩形回路（高度 $\to 0$），由 $\oint \mathbf{E} \cdot d\mathbf{l} = 0$：

$$E_{\text{above}}^{\parallel} \cdot l - E_{\text{below}}^{\parallel} \cdot l = 0$$

$$\boxed{E_{\text{above}}^{\parallel} = E_{\text{below}}^{\parallel}}$$

**电势的连续性**：

由 $V(\mathbf{b}) - V(\mathbf{a}) = -\int_\mathbf{a}^\mathbf{b} \mathbf{E} \cdot d\mathbf{l}$，当两点无限接近界面两侧时，积分趋于零：

$$\boxed{V_{\text{above}} = V_{\text{below}}}$$

电势在界面处连续（即使存在面电荷），但其**法向导数**不连续：

$$\frac{\partial V_{\text{above}}}{\partial n} - \frac{\partial V_{\text{below}}}{\partial n} = -\frac{\sigma}{\varepsilon_0}$$

---

## 2.4 静电学中的功与能

### 2.4.1 移动电荷的功

将试探电荷 $Q$ 在外电场 $\mathbf{E}$ 中从 $\mathbf{a}$ 搬运到 $\mathbf{b}$，外力需克服电场力做功：

$$W = -Q\int_\mathbf{a}^\mathbf{b} \mathbf{E} \cdot d\mathbf{l} = Q[V(\mathbf{b}) - V(\mathbf{a})]$$

因此，电势差就是单位电荷的势能差。电势的单位是伏特（V = J/C）。

### 2.4.2 点电荷系统的能量

**两个点电荷**：

将 $q_1$ 从无穷远带到位置 $\mathbf{r}_1$：不做功（空间中没有电场）。
将 $q_2$ 从无穷远带到位置 $\mathbf{r}_2$：需克服 $q_1$ 的电场做功：

$$W_2 = q_2 V_1(\mathbf{r}_2) = \frac{1}{4\pi\varepsilon_0}\frac{q_1 q_2}{r_{12}}$$

其中 $r_{12} = |\mathbf{r}_1 - \mathbf{r}_2|$。

**$N$ 个点电荷**：

逐个带入，第 $i$ 个电荷到来时要克服已有的所有电荷的电场：

$$W = \frac{1}{4\pi\varepsilon_0}\sum_{i=1}^{N}\sum_{j<i} \frac{q_i q_j}{r_{ij}}$$

利用对称性（每对计数两次，除以2）：

$$\boxed{W = \frac{1}{8\pi\varepsilon_0}\sum_{i=1}^{N}\sum_{j \neq i} \frac{q_i q_j}{r_{ij}} = \frac{1}{2}\sum_{i=1}^{N} q_i V(\mathbf{r}_i)}$$

其中 $V(\mathbf{r}_i)$ 是除 $q_i$ 以外的所有电荷在 $\mathbf{r}_i$ 处产生的电势。

**注意**：这个公式**不包含**电荷的自能（将点电荷从无穷远"组装"到一点所需的能量，对点电荷为无穷大）。

### 2.4.3 连续电荷分布的能量

将离散求和变为积分：

$$W = \frac{1}{2}\int \rho V \, d\tau$$

利用高斯定律 $\rho = \varepsilon_0 \nabla \cdot \mathbf{E}$ 和恒等式 $\nabla \cdot (V\mathbf{E}) = V(\nabla \cdot \mathbf{E}) + \mathbf{E} \cdot (\nabla V)$：

$$W = \frac{\varepsilon_0}{2}\int V(\nabla \cdot \mathbf{E}) \, d\tau = \frac{\varepsilon_0}{2}\int [\nabla \cdot (V\mathbf{E}) - \mathbf{E} \cdot \nabla V] \, d\tau$$

由 $\nabla V = -\mathbf{E}$，且 $\int \nabla \cdot (V\mathbf{E}) \, d\tau = \oint V\mathbf{E} \cdot d\mathbf{a} \to 0$（将积分推到无穷远，$V \sim 1/r$，$E \sim 1/r^2$，面积 $\sim r^2$，整体 $\sim 1/r \to 0$）：

$$\boxed{W = \frac{\varepsilon_0}{2}\int_{\text{全空间}} E^2 \, d\tau}$$

### 2.4.4 能量储存在哪里？——一个深刻的物理问题

我们现在有两个等价的能量表达式：

1. $W = \frac{1}{2}\int \rho V \, d\tau$ —— 能量在**电荷**所在处
2. $W = \frac{\varepsilon_0}{2}\int E^2 \, d\tau$ —— 能量在**电场**弥漫的全空间

**在静电学中**，这两个表达式给出相同的结果，无法区分。

**但在更一般的电动力学中**，我们会发现电磁波——不含任何电荷的纯电磁场——却携带着确定的能量和动量。这有力地支持了第二种观点：**能量储存在电场中**，能量密度为：

$$\boxed{u = \frac{\varepsilon_0}{2}E^2}$$

**细微差别**：公式 (1) 不包含自能（对点电荷有限），公式 (2) 包含自能（对点电荷发散）。对于连续分布，两者一致。点电荷自能的发散是经典电磁学的一个深层困难，直到量子电动力学（重整化）才获得解决。

---

**例题 2.10**：计算均匀带电球壳（半径 $R$，总电荷 $Q$）的静电能。

方法一（$W = \frac{\varepsilon_0}{2}\int E^2 d\tau$）：球壳内 $E = 0$，球壳外 $E = Q/(4\pi\varepsilon_0 r^2)$。

$$W = \frac{\varepsilon_0}{2}\int_R^\infty \left(\frac{Q}{4\pi\varepsilon_0 r^2}\right)^2 4\pi r^2 dr = \frac{Q^2}{8\pi\varepsilon_0}\int_R^\infty \frac{dr}{r^2} = \frac{Q^2}{8\pi\varepsilon_0 R}$$

方法二（$W = \frac{1}{2}\int \sigma V \, da$）：球面上 $V = Q/(4\pi\varepsilon_0 R)$，$\sigma = Q/(4\pi R^2)$。

$$W = \frac{1}{2}\cdot\frac{Q}{4\pi R^2}\cdot\frac{Q}{4\pi\varepsilon_0 R}\cdot 4\pi R^2 = \frac{Q^2}{8\pi\varepsilon_0 R}$$

两种方法给出相同结果。

$$\boxed{W = \frac{Q^2}{8\pi\varepsilon_0 R}}$$

---

## 2.5 导体

### 2.5.1 导体的基本性质

**导体**是含有大量自由电荷（通常是电子）的材料。在外电场的作用下，自由电荷迅速重新分布直到达到静电平衡。平衡态具有以下性质：

**性质 (i)：导体内部电场为零**

$$\mathbf{E} = 0 \quad \text{（导体内部）}$$

若不为零，自由电荷会受力运动，直到它们的重新分布消除内部电场——这就是"平衡"的含义。

**性质 (ii)：导体内部无净电荷**

由高斯定律 $\nabla \cdot \mathbf{E} = \rho/\varepsilon_0$，$\mathbf{E} = 0$ 立即给出 $\rho = 0$。

任何净电荷都必须分布在导体的**表面**上。

**性质 (iii)：导体表面是等势面**

由 $\mathbf{E} = -\nabla V = 0$（内部），$V$ 在导体内部为常数。由 $V$ 的连续性，表面上 $V$ 也取相同值。

等价地：沿导体表面 $\mathbf{E} \cdot d\mathbf{l} = 0$（因为 $\mathbf{E}$ 的切向分量为零），所以表面各点电势相等。

**性质 (iv)：导体表面的电场垂直于表面**

若 $\mathbf{E}$ 有切向分量，表面上的自由电荷会沿表面移动直到消除它。因此平衡时：

$$\mathbf{E} = \frac{\sigma}{\varepsilon_0}\hat{\mathbf{n}} \quad \text{（紧贴导体表面外侧）}$$

### 2.5.2 诱导电荷

当外电场施加于导体时，导体内的自由电荷重新分布，在表面形成**诱导电荷（Induced Charges）**。这些诱导电荷产生的电场恰好抵消外场在导体内部的效果。

**空腔内的场（法拉第笼效应）**：

考虑一个中空导体壳（内含空腔）。若空腔内没有电荷：

1. 导体内 $\mathbf{E} = 0$。
2. 取完全在导体内部的高斯面（包围空腔），$\oint \mathbf{E} \cdot d\mathbf{a} = 0$，故空腔内壁的总电荷为零。
3. 但仅仅总电荷为零还不够——内壁上**每一点**的面电荷密度都为零。

证明：假设空腔内壁存在面电荷。由于空腔内无自由电荷，$\oint \mathbf{E} \cdot d\mathbf{l} = 0$ 在空腔中成立，所以空腔中可以定义电势。空腔内的 $V$ 满足拉普拉斯方程 $\nabla^2 V = 0$，且在空腔壁上 $V = \text{const}$（因为壁是导体表面）。由拉普拉斯方程的极值性质（解不能在内部取极值），$V$ 在整个空腔中为常数，因此 $\mathbf{E} = -\nabla V = 0$。由边界条件 $\sigma = \varepsilon_0 E_\perp = 0$，内壁无电荷。

**物理后果**：导体壳完美屏蔽外部电场对内部的影响——这就是**法拉第笼**的原理。但注意：反过来不成立——腔内的电荷**会**在外表面诱导电荷，影响外部。

---

**例题 2.11**：一个带电荷 $+q$ 的点电荷位于接地导体球壳（半径 $R$）中心。求球壳内外的电场和球壳上的感应电荷。

**解**：由球对称，利用高斯定律。

- 球壳内（$r < R$）：以 $r < R$ 的球面为高斯面，$Q_{\text{enc}} = q$，故 $\mathbf{E} = \frac{q}{4\pi\varepsilon_0 r^2}\hat{\mathbf{r}}$。
- 球壳外（$r > R$）：球壳接地意味着 $V = 0$，等效于总电荷为零。内壁诱导 $-q$，外壁诱导 $+q$。但接地后外壁的 $+q$ 被"排走"（流入大地），外壁电荷为零。因此 $Q_{\text{enc}} = q + (-q) = 0$，$\mathbf{E} = 0$。

球壳内壁面电荷密度：$\sigma_{\text{inner}} = -q/(4\pi R^2)$。

### 2.5.3 导体表面的受力（静电压力）

导体表面的电荷处于电场中，会受到力。但这里有一个微妙之处：表面电荷"看到"的电场是什么？

表面正外方的电场是 $\sigma/\varepsilon_0$，正内方为 $0$。表面电荷本身感受到的场应取**上下两侧的平均值**：

$$\mathbf{E}_{\text{avg}} = \frac{1}{2}\left(\frac{\sigma}{\varepsilon_0} + 0\right)\hat{\mathbf{n}} = \frac{\sigma}{2\varepsilon_0}\hat{\mathbf{n}}$$

因此，单位面积上的力（**静电压力**）为：

$$\mathbf{f} = \sigma \mathbf{E}_{\text{avg}} = \frac{\sigma^2}{2\varepsilon_0}\hat{\mathbf{n}}$$

用电场强度表示（$\sigma = \varepsilon_0 E$，其中 $E$ 是紧贴表面外的场强）：

$$\boxed{f = \frac{\varepsilon_0}{2}E^2 = \frac{\sigma^2}{2\varepsilon_0}}$$

静电压力总是**向外**的——不论 $\sigma$ 的正负（因为 $\sigma^2 > 0$），电荷总是被电场向外拉。

**物理图景**：这个压力就是电场的能量密度 $u = \varepsilon_0 E^2/2$。这不是巧合——从能量角度看，如果导体表面外移 $dx$，释放出体积 $A \cdot dx$ 中储存的场能 $u \cdot A \cdot dx$，这恰好等于压力 $f$ 做的功 $f \cdot A \cdot dx$。

### 2.5.4 电容

**电容器**是一对导体，一个带电 $+Q$，另一个带电 $-Q$。两导体之间的电势差 $V = V_+ - V_-$ 与 $Q$ 成正比（因为 $\mathbf{E}$ 正比于 $Q$，$V$ 正比于 $\mathbf{E}$）：

$$\boxed{C = \frac{Q}{V}}$$

$C$ 称为**电容**，单位为法拉（F = C/V）。电容只取决于导体的几何形状，与 $Q$ 和 $V$ 无关。

---

**例题 2.12**：平行板电容器。

两块面积为 $A$、间距为 $d$ 的平行导体板。板间电场（忽略边缘效应）：

$$E = \frac{\sigma}{\varepsilon_0} = \frac{Q}{\varepsilon_0 A}$$

电势差：$V = Ed = Qd/(\varepsilon_0 A)$。

$$\boxed{C = \frac{\varepsilon_0 A}{d}}$$

**物理直觉**：面积越大，可容纳的电荷越多；间距越小，同样电荷产生的电势差越小——两者都增大电容。

---

**例题 2.13**：同轴圆柱电容器。

内圆柱半径 $a$，外圆柱半径 $b$，长度 $L$（$L \gg b$）。两柱之间电场（利用高斯定律）：

$$\mathbf{E} = \frac{\lambda}{2\pi\varepsilon_0 s}\hat{\mathbf{s}} = \frac{Q}{2\pi\varepsilon_0 L s}\hat{\mathbf{s}}$$

电势差：
$$V = -\int_b^a \mathbf{E} \cdot d\mathbf{l} = -\int_b^a \frac{Q}{2\pi\varepsilon_0 L s}ds = \frac{Q}{2\pi\varepsilon_0 L}\ln\frac{b}{a}$$

$$\boxed{C = \frac{2\pi\varepsilon_0 L}{\ln(b/a)}}$$

---

**例题 2.14**：同心球壳电容器。

内球壳半径 $a$，外球壳半径 $b$。

$$V = -\int_b^a \frac{Q}{4\pi\varepsilon_0 r^2}dr = \frac{Q}{4\pi\varepsilon_0}\left(\frac{1}{a} - \frac{1}{b}\right)$$

$$\boxed{C = 4\pi\varepsilon_0\frac{ab}{b-a}}$$

当 $b \to \infty$ 时，$C = 4\pi\varepsilon_0 a$——这是孤立导体球的电容（以无穷远为"另一极板"）。

---

**电容器储存的能量**：

$$W = \frac{1}{2}CV^2 = \frac{Q^2}{2C} = \frac{1}{2}QV$$

---

## 习题

### 基础计算

**2.1** 两个电荷 $+q$ 和 $-q$ 位于 $x$ 轴上，分别在 $x = +d/2$ 和 $x = -d/2$。
(a) 求 $z$ 轴上任意点 $(0, 0, z)$ 处的电场（精确结果）。
(b) 对 $z \gg d$ 的远区取近似，验证结果与偶极子公式一致。
(c) 求 $x$ 轴上（$x > d/2$）的电场。

**2.2** 一根半径为 $R$ 的无限长圆柱体携带体电荷密度 $\rho = ks$（$s$ 是到轴线的距离，$k$ 为常数）。利用高斯定律求柱内外的电场。

**2.3** 一个厚度为 $2d$ 的无限大平板携带均匀体电荷密度 $\rho$。以平板中心为原点，求电场 $\mathbf{E}$ 作为 $y$ 坐标的函数，并画出 $E$ 对 $y$ 的图。

**2.4** 用两种方法计算均匀带电实心球（半径 $R$，总电荷 $Q$）的静电能量 $W$：
(a) 用 $W = \frac{\varepsilon_0}{2}\int E^2 d\tau$。
(b) 用"组装"方法：逐层从无穷远搬运球壳电荷。验证两种方法结果一致。

**2.5** 一个半径为 $R$ 的导体球携带总电荷 $Q$。求"北半球"与"南半球"之间的斥力。（提示：利用静电压力。）

### 概念理解

**2.6** **高斯定律的局限性**：为什么高斯定律对任何电荷分布都成立，但作为计算工具只在高对称性问题中有用？用一个有限长带电线段为例说明：写出高斯定律，但解释为什么无法直接求解 $\mathbf{E}$。

**2.7** **电势的相对性**：电势的零点可以任意选取（如同海拔高度的零点）。解释为什么**电势差** $V(\mathbf{b}) - V(\mathbf{a})$ 是物理上有意义的量，而绝对电势本身不是。在什么条件下，我们习惯取无穷远为零势点？什么时候这个选择不方便（提示：无限长带电线）？

**2.8** **壳层定理的深层原因**：利用高斯定律证明：均匀带电球壳内部电场为零。然后解释：如果库仑定律不是精确的平方反比（如 $1/r^{2+\epsilon}$），球壳内是否仍有 $\mathbf{E} = 0$？这与高斯定律的什么性质有关？

**2.9** **能量悖论**：考虑一个点电荷 $q$。利用 $W = \frac{\varepsilon_0}{2}\int E^2 d\tau$ 计算其能量，说明结果为什么发散。这对经典电动力学意味着什么？（提示：思考"点粒子"假设的合理性。）

### 拓展应用

**2.10** **叠加法的威力**：一个球形空腔被挖去了一个偏心球体（球心不重合）。设大球半径 $R$，体电荷密度 $\rho$，球形空腔半径 $a$，球心偏移 $\mathbf{d}$。证明空腔内电场均匀，且 $\mathbf{E} = \frac{\rho}{3\varepsilon_0}\mathbf{d}$。

**提示**：将此问题视为两个均匀带电球的叠加——一个正电荷球（大球）和一个负电荷球（空腔）。

**2.11** **编程练习**：编写Python程序，利用库仑定律数值计算并可视化以下电荷分布的电场线和等势线：
(a) 一个点电荷。
(b) 一个电偶极子（两个等大反号点电荷）。
(c) 四个点电荷组成的四极子（$+q, -q, +q, -q$ 在正方形顶点）。

**提示**：在二维网格上计算 $\mathbf{E}$ 的分量，使用 `matplotlib.pyplot.streamplot` 画场线，`contour` 画等势线。

```python
import numpy as np
import matplotlib.pyplot as plt

# 示例框架：点电荷的电场
def E_point_charge(q, r0, x, y):
    """计算点电荷 q 在位置 r0=(x0,y0) 产生的电场"""
    dx = x - r0[0]
    dy = y - r0[1]
    r = np.sqrt(dx**2 + dy**2)
    r = np.maximum(r, 1e-10)  # 避免除零
    Ex = q * dx / r**3
    Ey = q * dy / r**3
    return Ex, Ey

# 创建网格
x = np.linspace(-3, 3, 200)
y = np.linspace(-3, 3, 200)
X, Y = np.meshgrid(x, y)

# 偶极子: +q 在 (0, 0.5), -q 在 (0, -0.5)
Ex1, Ey1 = E_point_charge(1, (0, 0.5), X, Y)
Ex2, Ey2 = E_point_charge(-1, (0, -0.5), X, Y)
Ex, Ey = Ex1 + Ex2, Ey1 + Ey2

# 限制场强避免奇点附近的图形溢出
E_mag = np.sqrt(Ex**2 + Ey**2)
Ex_norm = Ex / E_mag
Ey_norm = Ey / E_mag

fig, ax = plt.subplots(figsize=(8, 8))
ax.streamplot(X, Y, Ex_norm, Ey_norm, density=2, color=np.log(E_mag), cmap='inferno')
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_aspect('equal')
ax.set_title('Electric Dipole Field Lines')
plt.savefig('dipole_field_lines.png', dpi=150)
plt.show()
```

### Griffiths 教材精选习题

> 以下习题直接对应 Griffiths《Introduction to Electrodynamics》第五版题号，标注为 **[G x.xx]**。

**[G 2.12] 球壳电场——高斯定律的经典应用**

利用高斯定律求均匀带电球壳（半径 $R$，面电荷密度 $\sigma$）内外的电场。与直接积分法（对球面上的每一面元计算其贡献再叠加）相比，体会高斯定律的高效。

---

**[G 2.22] 均匀带电实心球的电势分布**

求均匀带电实心球（半径 $R$，总电荷 $q$）内外的电势 $V(r)$。以无穷远为参考点。分别计算各区域 $V$ 的梯度，验证其给出正确的电场。画出 $V(r)$ 的示意图。

---

**[G 2.31] 电场边界条件验证**

(a) 验证例题 2.5、例题 2.6 和 Prob 2.12 的结果与边界条件（$E_\perp$ 的跃变 $= \sigma/\varepsilon_0$）一致。

(b) 利用高斯定律求无限长中空圆柱管（面电荷密度 $\sigma$）内外的电场，验证其满足边界条件。

(c) 验证均匀带电球壳（例题 2.8）的结果与 $E_\perp$ 和 $E_\parallel$ 的边界条件一致。

---

**[G 2.42] 平行板电容器的静电压力**

两块面积为 $A$ 的大金属板相距 $d$，每块板上放置电荷 $Q$。求板上的静电压力。

**提示**：注意两板携带同号电荷，与标准电容器（$+Q$ 和 $-Q$）不同。先求板间和板外的电场（利用叠加），然后用 $f = \sigma^2/(2\varepsilon_0)$ 或直接从场的能量密度出发计算压力。

---

**[G 2.44] 同轴圆柱电缆的电容**

求两个同轴金属圆柱管（半径分别为 $a$ 和 $b$，$a < b$）单位长度的电容。这是同轴电缆的基本参数，其值决定了电缆的信号传输特性。

---

**[G 2.52] 均匀带电圆盘边缘的电势——精巧的积分**

求均匀带电圆盘（半径 $R$，面电荷密度 $\sigma$）边缘处的电势。

**提示**：首先用量纲分析说明 $V = k\sigma R/(\pi\varepsilon_0)$，其中 $k$ 是某个无量纲常数。将 $k$ 表示为一个定积分，然后尝试解析计算（若能）或用数值方法求解。

---

## Key Takeaway（本章核心要点）

1. **场是物理实在，而非计算工具**：电场 $\mathbf{E}(\mathbf{r})$ 在空间每一点都有确定的值，携带能量和动量。它不仅仅是计算力的中间步骤——它是物质世界的一部分。

2. **高斯定律连接了源与场**：$\nabla \cdot \mathbf{E} = \rho/\varepsilon_0$ 将电场的散度与电荷密度直接挂钩。对称性问题中，它是最有力的计算工具。但请记住：高斯定律只给出了 $\mathbf{E}$ 的散度；要完全确定 $\mathbf{E}$，还需要旋度条件 $\nabla \times \mathbf{E} = 0$（在静电学中）。

3. **电势将矢量问题化为标量问题**：$\mathbf{E} = -\nabla V$ 使我们可以先计算标量 $V$（更简单的积分），再求梯度得到 $\mathbf{E}$。泊松方程 $\nabla^2 V = -\rho/\varepsilon_0$ 是静电学的核心微分方程。

4. **边界条件是求解问题的钥匙**：电场法向分量在面电荷处跃变 $\sigma/\varepsilon_0$，切向分量连续，电势连续。这些条件在第3章的边值问题中将发挥决定性作用。

5. **能量储存在场中**：$u = \varepsilon_0 E^2/2$ 是电场的能量密度。这一观点在第9章（电磁波携带能量）和第8章（坡印廷定理）中将得到最充分的验证。

6. **导体是"被电场塑造的物体"**：导体内部 $\mathbf{E} = 0$、电荷分布在表面、表面等势——这些性质不是假设，而是自由电荷对电场响应的**自然结果**。导体的行为完全由静电平衡决定。

---

**致读者**：本章的方法论可以概括为一个策略升级路径：**库仑定律（通用但笨拙）→ 高斯定律（优雅但需要对称性）→ 电势（万能但需要积分/求解微分方程）**。在第3章中，你将学到更强大的技巧——镜像法、分离变量法和多极展开——来攻克没有简单对称性的边值问题。那时你会发现，本章建立的电势理论和边界条件是一切技巧的合法性基础。
