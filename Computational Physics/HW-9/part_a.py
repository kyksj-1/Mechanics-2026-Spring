"""
Part A: 2x2 矩阵对角化

H = [[a, c], [c*, -a]]

解析推导 + 数值验证 e^{iHt} 的闭式表达.
"""

import numpy as np


def hamiltonian_2x2(a: float, c: complex) -> np.ndarray:
    return np.array([[a, c], [np.conj(c), -a]], dtype=np.complex128)


def analytic_exp_iHt(a: float, c: complex, t: float) -> np.ndarray:
    """e^{iHt} 的解析闭式: cos(E t) I + i sin(E t) / E * H,  E = sqrt(a^2 + |c|^2)."""
    E = np.sqrt(a ** 2 + np.abs(c) ** 2)
    H = hamiltonian_2x2(a, c)
    I = np.eye(2, dtype=np.complex128)
    if E == 0:
        return I
    return np.cos(E * t) * I + 1j * np.sin(E * t) / E * H


def numeric_exp_iHt(a: float, c: complex, t: float) -> np.ndarray:
    """对 H 数值对角化, 用 Q e^{iLambda t} Q^{-1} 计算."""
    H = hamiltonian_2x2(a, c)
    evals, Q = np.linalg.eigh(H)  # 厄米 -> Q^{-1} = Q^H
    Lambda = np.diag(np.exp(1j * evals * t))
    return Q @ Lambda @ Q.conj().T


def main():
    print("=" * 70)
    print("Part A: 二能级哈密顿量 H = [[a, c], [c*, -a]] 的对角化与时间演化")
    print("=" * 70)

    # 取一组示例参数验证: a = 0.7, c = 0.3 + 0.4j
    a = 0.7
    c = 0.3 + 0.4j
    H = hamiltonian_2x2(a, c)
    print(f"\n示例参数: a = {a}, c = {c}")
    print(f"|c| = {np.abs(c):.6f}")
    print(f"E = sqrt(a^2 + |c|^2) = {np.sqrt(a ** 2 + np.abs(c) ** 2):.6f}")

    # 数值对角化
    evals, Q = np.linalg.eigh(H)
    print(f"\n特征值 Lambda = diag({evals[0]:.6f}, {evals[1]:.6f})")
    print(f"理论值 ±E   = ±{np.sqrt(a ** 2 + np.abs(c) ** 2):.6f}")

    print(f"\n本征向量矩阵 Q =")
    print(Q)

    # 验证 Q^{-1} = Q^H
    Qinv = np.linalg.inv(Q)
    QH = Q.conj().T
    print(f"\n|| Q^{{-1}} - Q^H ||_F = {np.linalg.norm(Qinv - QH):.3e}")
    print("Q 是酉矩阵 (H 厄米 -> 本征向量可取正交归一 -> Q^{-1} = Q^H)")

    # 验证 H = Q Lambda Q^H
    H_reconstructed = Q @ np.diag(evals) @ QH
    print(f"\n|| H - Q Lambda Q^H ||_F = {np.linalg.norm(H - H_reconstructed):.3e}")

    # 验证 e^{iHt} 解析公式
    print("\n--- 验证 e^{iHt} 的闭式公式 ---")
    print("公式: e^{iHt} = cos(Et) I + i sin(Et)/E * H,   E = sqrt(a^2 + |c|^2)")
    for t in [0.5, 1.3, 3.7, 10.0]:
        U_ana = analytic_exp_iHt(a, c, t)
        U_num = numeric_exp_iHt(a, c, t)
        err = np.linalg.norm(U_ana - U_num)
        unit_err = np.linalg.norm(U_ana.conj().T @ U_ana - np.eye(2))
        print(f"  t = {t:6.2f}:  || ana - num ||_F = {err:.3e},  unitarity err = {unit_err:.3e}")

    # 极端情况: c = 0
    print("\n--- 退化情形 c = 0 ---")
    print("此时 H 已是对角的, e^{iHt} = diag(e^{iat}, e^{-iat})")
    a2, c2, t2 = 1.5, 0.0 + 0.0j, 2.0
    U_ana = analytic_exp_iHt(a2, c2, t2)
    U_diag = np.diag([np.exp(1j * a2 * t2), np.exp(-1j * a2 * t2)])
    print(f"  || ana - diag(e^{{iat}}, e^{{-iat}}) || = {np.linalg.norm(U_ana - U_diag):.3e}")

    print("\nPart A 验证完成.")


if __name__ == "__main__":
    main()
