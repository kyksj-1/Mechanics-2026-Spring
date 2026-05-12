# A2. MCMC 细致平衡与 Metropolis 设计

## 一般细致平衡
对目标平衡分布 \(\pi(C)\)（构型记作 \(C\)），马尔可夫链转移概率写作
\[
W(C\to C') = g(C\to C')A(C\to C')
\]
其中 \(g\) 是提议概率，\(A\) 是接受概率。细致平衡要求
\[
\pi(C)W(C\to C')=\pi(C')W(C'\to C).
\]
因此可取
\[
\frac{A(C\to C')}{A(C'\to C)}=
\frac{\pi(C')g(C'\to C)}{\pi(C)g(C\to C')}.
\]

## Ising 模型中的权重
本题中
\[
H(C)=-\sum_{\langle ij\rangle}\sigma_i\sigma_j,
\quad
\pi(C)=\frac{e^{-\beta H(C)}}{Z}.
\]
若一次更新只尝试翻转一个格点，构型变化为 \(C\to C'\)，能量差
\[
\Delta E=H(C')-H(C)=2\sigma_i\sum_{j\in nn(i)}\sigma_j.
\]

## 过程与逆过程
- 过程：在 \(L^2\) 个格点中等概率随机选一个点 \(i\)，尝试翻转 \(\sigma_i\to-\sigma_i\)。
- 逆过程：在新构型中再次选择同一个点 \(i\) 翻回去。

因为选点是均匀随机，
\[
g(C\to C')=g(C'\to C)=\frac{1}{L^2}.
\]
所以
\[
\frac{A(C\to C')}{A(C'\to C)}=e^{-\beta\Delta E}.
\]
采用 Metropolis 选择
\[
A(C\to C')=\min\{1,e^{-\beta\Delta E}\}.
\]
这保证细致平衡成立。
