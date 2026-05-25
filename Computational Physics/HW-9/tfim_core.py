"""
横场 Ising 模型 (Transverse-field Ising model, TFIM) 核心工具模块

H = -J * sum_i sigma_i^z sigma_{i+1}^z - h * sum_i sigma_i^x

约定:
    基矢使用计算基 |s_{L-1} ... s_1 s_0>, 其中比特 i=0 对应整数最低位.
    sigma_i^z |s> = (1 - 2*s_i) |s>     (s_i=0 -> +1; s_i=1 -> -1)
    sigma_i^x |s> = |s XOR (1<<i)>

实现亮点:
    - 全部用 numpy 位运算向量化构造稀疏 H, 避免 Python 级 for 循环.
    - 对角部分用对角矩阵, 非对角部分用 CSR 一次性构建.
    - L=18 (dim=262144) 仍可在内存中构造并参与 Krylov 演化.
"""

import numpy as np
from scipy import sparse


def build_tfim_hamiltonian(L: int, J: float, h: float, pbc: bool = True) -> sparse.csr_matrix:
    """构造 L 格点 TFIM 的稀疏哈密顿量.

    参数
    ----
    L : 格点数
    J, h : 哈密顿量参数
    pbc : 是否使用周期边界条件
    """
    N = 1 << L  # 2**L
    states = np.arange(N, dtype=np.int64)

    # 对角的 ZZ 部分: -J * sum sigma_z_i sigma_z_{i+1}
    # sigma_z_i 在基态 |s> 上的本征值为 (1-2*((s>>i)&1))
    diag = np.zeros(N, dtype=np.float64)
    pairs = [(i, (i + 1) % L) for i in range(L)] if pbc else [(i, i + 1) for i in range(L - 1)]
    for i, j in pairs:
        zi = 1 - 2 * ((states >> i) & 1)
        zj = 1 - 2 * ((states >> j) & 1)
        diag += -J * zi * zj

    # 非对角的 X 部分: -h * sum sigma_x_i, 把比特 i 翻转, 矩阵元 = -h
    # 总非零数 = L * N, 每行 L 个
    if h != 0.0:
        rows = np.empty(L * N, dtype=np.int64)
        cols = np.empty(L * N, dtype=np.int64)
        for i in range(L):
            flipped = states ^ (1 << i)
            rows[i * N:(i + 1) * N] = states
            cols[i * N:(i + 1) * N] = flipped
        data = np.full(L * N, -h, dtype=np.float64)
        H_off = sparse.csr_matrix((data, (rows, cols)), shape=(N, N))
    else:
        H_off = sparse.csr_matrix((N, N), dtype=np.float64)

    H_diag = sparse.diags(diag, format='csr')
    return (H_diag + H_off).tocsr()


def sigma_x_op(L: int, site: int) -> sparse.csr_matrix:
    """构造 sigma_x 算符在第 site 个格点上 (site=0 为最低位)."""
    N = 1 << L
    states = np.arange(N, dtype=np.int64)
    flipped = states ^ (1 << site)
    data = np.ones(N, dtype=np.float64)
    return sparse.csr_matrix((data, (states, flipped)), shape=(N, N))


def sigma_z_op(L: int, site: int) -> sparse.dia_matrix:
    """构造 sigma_z 算符在第 site 个格点上, 对角矩阵."""
    N = 1 << L
    states = np.arange(N, dtype=np.int64)
    diag = (1 - 2 * ((states >> site) & 1)).astype(np.float64)
    return sparse.diags(diag, format='dia')


def expectation(psi: np.ndarray, op: sparse.spmatrix) -> complex:
    """计算 <psi| op |psi>, psi 为列向量."""
    return np.vdot(psi, op @ psi)


def neel_state_index(L: int) -> int:
    """Néel 态 |↑↓↑↓...> 的整数下标.

    约定: 偶数 site (0,2,...) = up (0), 奇数 site (1,3,...) = down (1).
    这对应 "奇数格点朝上, 偶数格点朝下" 的 1-indexed 表述
    (题面 site 从 1 起算, 奇数 site 1,3,... 朝上;
     代码内 0-indexed 即偶数 site 0,2,... 朝上).
    """
    idx = 0
    for i in range(L):
        if i % 2 == 1:  # 奇数 site (0-indexed) -> 题面偶数 site -> 朝下 -> 比特=1
            idx |= (1 << i)
    return idx


def neel_state_vector(L: int) -> np.ndarray:
    """Néel 初态的列向量 (在计算基中)."""
    N = 1 << L
    psi = np.zeros(N, dtype=np.complex128)
    psi[neel_state_index(L)] = 1.0
    return psi
