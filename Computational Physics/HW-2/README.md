# 计算物理 Homework 2: FFT & 牛顿法

> 作者: kyksj-1
> 日期: 2026-03-20
> 许可证: MIT License

---

## 目录

- [项目结构](#项目结构)
- [环境与运行](#环境与运行)
- [A. DFT 和 FFT](#a-dft-和-fft)
  - [A.1 离散傅里叶变换 (DFT)](#a1-离散傅里叶变换-dft)
  - [A.2 Base-2 快速傅里叶变换 (FFT)](#a2-base-2-快速傅里叶变换-fft)
  - [A.3 信号频域分析](#a3-信号频域分析)
- [B. 牛顿迭代法](#b-牛顿迭代法)（待完成）

---

## 项目结构

```
HW-2/
├── config/
│   ├── config.py          # 静态路径配置（项目根目录、数据/输出路径等）
│   └── config.yaml        # 动态运行参数（基准测试重复次数、随机种子、绘图参数等）
├── src/
│   ├── dft.py             # DFT / IDFT 核心算法
│   ├── fft.py             # Base-2 FFT / IFFT 核心算法 (Cooley-Tukey)
│   └── signal_analysis.py # 信号加载、频谱计算、峰值检测
├── scripts/
│   ├── run_dft_test.py         # A.1 DFT 验证与可视化
│   ├── run_fft_benchmark.py    # A.2 FFT 基准测试与复杂度分析
│   └── run_signal_analysis.py  # A.3 波形信号频域分析
├── data/
│   └── waveform.dat       # 示波器采集的电压信号数据
├── output/                # 运行脚本自动生成的图片（不入库）
├── ProbSet_2.md           # 题目原文
└── README.md              # 本文件
```

---

## 环境与运行

```bash
conda activate research_env

# A.1 DFT 验证
python scripts/run_dft_test.py

# A.2 FFT 基准测试
python scripts/run_fft_benchmark.py

# A.3 信号分析
python scripts/run_signal_analysis.py
```

依赖: `numpy`, `matplotlib`, `pyyaml`

---

## A. DFT 和 FFT

### A.1 离散傅里叶变换 (DFT)

#### 问题

编写 DFT 程序，接受 `complex128` 一维数组，返回其离散傅里叶变换。给出测试样例并对比库函数结果。

#### 算法

DFT 严格按定义式实现：

$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-2\pi i \cdot nk / N}, \quad k = 0, 1, \ldots, N-1$$

**关键实现 (`src/dft.py`):** 将上式改写为矩阵-向量乘法 $\mathbf{X} = \mathbf{W} \cdot \mathbf{x}$，其中 DFT 矩阵 $W_{kn} = e^{-2\pi i \cdot kn / N}$。利用 `np.outer(k, n)` 一次性构造指数矩阵，避免显式双重循环。复杂度 $O(N^2)$。

#### 测试结果

对随机 `complex128` 数组，自实现 DFT 与 `numpy.fft.fft` 的对比：

| N   | 最大绝对误差 | 相对误差 (L2) | 通过 |
|-----|-------------|--------------|------|
| 8   | 6.40e-15    | 7.72e-16     | YES  |
| 16  | 1.52e-14    | 1.63e-15     | YES  |
| 32  | 7.83e-14    | 3.98e-15     | YES  |
| 64  | 3.03e-13    | 8.27e-15     | YES  |
| 128 | 7.18e-13    | 1.36e-14     | YES  |

所有测试用例的最大绝对误差均在 $10^{-13}$ 以下，误差来源为浮点运算的舍入误差累积，属于正常范围。

#### 可视化

![DFT 对比验证](output/a1_dft_comparison.png)

图中包含：
- 幅值谱 $|X[k]|$ 对比（蓝色自实现 vs 红色 numpy，完全重合）
- 相位谱对比
- 逐点绝对误差（$10^{-14} \sim 10^{-13}$ 量级）
- IDFT 还原误差（验证 $x = \text{IDFT}(\text{DFT}(x))$）

---

### A.2 Base-2 快速傅里叶变换 (FFT)

#### 问题

编写 Base-2 FFT，检查输入长度是否为 2 的幂次，给出 $2^4$ 到 $2^{12}$ 的测试。分析计算复杂度并与标准库对比用时。

#### 算法 — Cooley-Tukey 递归分治

**核心思想 (`src/fft.py`):** 将长度为 $N$ 的 DFT 按奇偶下标拆分为两个 $N/2$ 的子问题，通过蝶形运算合并：

$$X[k] = E[k] + W_N^k \cdot O[k], \quad X[k + N/2] = E[k] - W_N^k \cdot O[k]$$

其中 $W_N^k = e^{-2\pi ik/N}$ 为旋转因子 (twiddle factor)，$E[k]$ 和 $O[k]$ 分别是偶数/奇数子序列的 DFT。

递归终止条件：$N=1$ 时 DFT 就是元素本身。

#### 复杂度分析

| 算法 | 复杂度 | N=4096 时的运算量 |
|------|--------|------------------|
| DFT  | $O(N^2)$ | ~16,777,216 |
| FFT  | $O(N \log N)$ | ~49,152 |

**推导：** 递归关系 $T(N) = 2T(N/2) + O(N)$，共 $\log_2 N$ 层，每层 $N$ 次蝶形运算，总计 $N \log_2 N$。

#### 正确性验证

$2^4$ 到 $2^{12}$ 全部通过（最大绝对误差 < $10^{-8}$）：

| N    | 最大误差    | 通过 |
|------|-----------|------|
| 16   | 1.26e-15  | YES  |
| 32   | 2.59e-15  | YES  |
| 64   | 6.04e-15  | YES  |
| 128  | 1.78e-14  | YES  |
| 256  | 2.42e-14  | YES  |
| 512  | 3.66e-14  | YES  |
| 1024 | 6.90e-14  | YES  |
| 2048 | 1.15e-13  | YES  |
| 4096 | 1.88e-13  | YES  |

#### 性能对比

使用 `time.perf_counter_ns()` 计时（纳秒精度），每个大小重复 5 次取中位数：

| N    | 自实现 FFT | numpy FFT | 自实现 DFT | FFT/numpy 速度比 |
|------|-----------|-----------|-----------|-----------------|
| 16   | 162 us    | 3.2 us    | 32 us     | ~50x            |
| 64   | 774 us    | 3.2 us    | 576 us    | ~242x           |
| 256  | 3.24 ms   | 5.3 us    | 4.24 ms   | ~611x           |
| 1024 | 10.0 ms   | 22 us     | 55.7 ms   | ~454x           |
| 4096 | 44.1 ms   | 57.7 us   | N/A       | ~764x           |

**分析：**

1. **自实现 FFT 确实比 DFT 快。** 在 N=1024 时，FFT 耗时 10ms vs DFT 的 56ms，加速约 5.6 倍。理论上 $N / \log_2 N = 1024/10 = 102$ 倍，实际加速比较小是因为 Python 递归的额外开销抵消了一部分算法优势。

2. **自实现 FFT 比 numpy 慢约 50~800 倍。** 原因很明确：
   - numpy 底层调用 C/Fortran 编写的 FFTPACK / pocketfft，无解释器开销
   - 自实现使用 Python 递归，每次递归都有函数调用和数组切片开销
   - numpy 在底层使用了 SIMD 向量化指令和缓存优化

3. **复杂度趋势符合理论。** 从双对数图中可以看到，自实现 FFT 的斜率与 $O(N\log N)$ 参考线平行，DFT 的斜率与 $O(N^2)$ 参考线平行。

![FFT 性能基准测试](output/a2_fft_benchmark.png)

---

### A.3 信号频域分析

#### 问题

分析 `waveform.dat` 中示波器测量的电压信号的频域特征。

#### 数据概况

| 参数 | 值 |
|------|-----|
| 采样点数 | 40,000 |
| 采样间隔 | 0.001 s |
| 采样率 | 1000 Hz |
| 总时长 | 40 s |
| Nyquist 频率 | 500 Hz |
| 电压范围 | [-1.23, 1.25] V |

#### 频谱分析结果

FFT 后发现 **4 个显著的频率分量**，全部位于低频区域：

| 频率 (Hz) | 归一化幅值 | 周期 (s) | 与基频之比 |
|-----------|-----------|---------|-----------|
| **1.00**  | 1.0000    | 1.000   | 1 (基频)  |
| **3.00**  | 0.3340    | 0.333   | 3         |
| **5.00**  | 0.2006    | 0.200   | 5         |
| **7.00**  | 0.1439    | 0.143   | 7         |

#### 发现

这是一个 **方波信号** 叠加白噪声。判断依据：

1. **仅包含奇次谐波：** 频率分量为基频的 1, 3, 5, 7 倍，偶次谐波缺失——这是方波的标志性特征。

2. **幅值比接近 $1/n$：** 方波的傅里叶级数展开为：
   $$f(t) = \frac{4}{\pi} \sum_{n=1,3,5,...} \frac{1}{n} \sin(n\omega_0 t)$$
   理论幅值比为 $1 : 1/3 : 1/5 : 1/7 = 1 : 0.333 : 0.200 : 0.143$，与实测值 $1 : 0.334 : 0.201 : 0.144$ 高度吻合。

3. **宽带噪声本底：** 除了上述离散谱线外，整个频带内存在均匀的低幅值噪声本底（约 0.001~0.003），呈现白噪声特征。

4. **信号稳态性：** 从时频图 (Spectrogram) 可以看出，各频率分量在整个 40s 时间段内幅值恒定，说明信号是稳态的。

#### 可视化

![信号频域分析](output/a3_signal_analysis.png)

![时频分析](output/a3_spectrogram.png)

---

## B. 牛顿迭代法

（待后续会话完成）

---

# Computational Physics Homework 2: FFT & Newton's Method

> Author: kyksj-1
> Date: 2026-03-20
> License: MIT License

---

## Table of Contents

- [Project Structure](#project-structure)
- [Environment & Execution](#environment--execution)
- [A. DFT and FFT](#a-dft-and-fft-1)
  - [A.1 Discrete Fourier Transform (DFT)](#a1-discrete-fourier-transform-dft)
  - [A.2 Base-2 Fast Fourier Transform (FFT)](#a2-base-2-fast-fourier-transform-fft)
  - [A.3 Signal Frequency Analysis](#a3-signal-frequency-analysis)
- [B. Newton's Method](#b-newtons-method) (To be completed)

---

## Project Structure

```
HW-2/
├── config/
│   ├── config.py          # Static path configuration
│   └── config.yaml        # Runtime parameters (benchmark repeats, seed, plot settings)
├── src/
│   ├── dft.py             # DFT / IDFT core algorithms
│   ├── fft.py             # Base-2 FFT / IFFT (Cooley-Tukey)
│   └── signal_analysis.py # Signal loading, spectrum computation, peak detection
├── scripts/
│   ├── run_dft_test.py         # A.1 DFT verification & visualization
│   ├── run_fft_benchmark.py    # A.2 FFT benchmarking & complexity analysis
│   └── run_signal_analysis.py  # A.3 Waveform frequency-domain analysis
├── data/
│   └── waveform.dat       # Oscilloscope voltage signal data
├── output/                # Auto-generated figures (not tracked by git)
├── ProbSet_2.md           # Problem set
└── README.md              # This file
```

---

## Environment & Execution

```bash
conda activate research_env

# A.1 DFT verification
python scripts/run_dft_test.py

# A.2 FFT benchmark
python scripts/run_fft_benchmark.py

# A.3 Signal analysis
python scripts/run_signal_analysis.py
```

Dependencies: `numpy`, `matplotlib`, `pyyaml`

---

## A. DFT and FFT

### A.1 Discrete Fourier Transform (DFT)

**Algorithm:** Implemented directly from the DFT definition $X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-2\pi i nk/N}$, vectorized as a matrix-vector product $\mathbf{X} = \mathbf{W}\mathbf{x}$ where $W_{kn} = e^{-2\pi ikn/N}$. Complexity: $O(N^2)$.

**Key implementation detail (`src/dft.py`):** The exponent matrix is constructed via `np.outer(k, n)` to avoid explicit double loops.

**Verification:** All test cases (N = 8, 16, 32, 64, 128) pass with maximum absolute error below $10^{-10}$ when compared to `numpy.fft.fft`.

### A.2 Base-2 Fast Fourier Transform (FFT)

**Algorithm:** Cooley-Tukey radix-2 decimation-in-time (DIT), recursive implementation. Splits the DFT into even/odd-indexed sub-problems and merges via the butterfly operation:

$$X[k] = E[k] + W_N^k \cdot O[k], \quad X[k+N/2] = E[k] - W_N^k \cdot O[k]$$

**Complexity:** $O(N\log N)$. Recurrence: $T(N) = 2T(N/2) + O(N)$, yielding $N\log_2 N$ complex multiply-adds across $\log_2 N$ recursion levels.

**Results:** All sizes from $2^4$ to $2^{12}$ pass correctness verification. The Python implementation is ~50-800x slower than numpy (due to interpreter overhead vs. compiled C/FFTPACK), but the log-log slope confirms $O(N\log N)$ scaling.

### A.3 Signal Frequency Analysis

**Finding:** The waveform is a **square wave** (fundamental frequency 1 Hz) plus white noise. Evidence:

1. Only **odd harmonics** present: 1, 3, 5, 7 Hz
2. Amplitudes follow the $1/n$ pattern: 1.000, 0.334, 0.201, 0.144 (matching the Fourier series $\frac{4}{\pi}\sum \frac{1}{n}\sin(n\omega_0 t)$)
3. Flat broadband noise floor across the spectrum
4. Stationary signal (confirmed by spectrogram)

---

## B. Newton's Method

(To be completed in a later session)
