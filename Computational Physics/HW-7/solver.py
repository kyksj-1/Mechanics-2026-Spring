import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

class FrustratedIsing2D:
    """
    A class to simulate and sample the ground states of the 2D frustrated Ising model.
    A sublattice has horizontal A-A bonds.
    B sublattice has vertical B-B bonds.
    A and B have 4 cross bonds (forming squares).
    """

    def __init__(self, L=16):
        """
        Initialize the Ising model on L x L cells.
        Each cell has 1 A site and 1 B site, total 2L^2 sites.
        """
        self.L = L
        if L % 2 != 0:
            print("Warning: L is recommended to be even for periodic ground states.")
            
    def generate_ground_state(self):
        """
        Generate a random ground state.
        Rule: 
        A chains (along x) are antiferromagnetic.
        B chains (along y) are antiferromagnetic.
        A spins only depend on y layer: A(x,y) = a_y * (-1)^x
        B spins only depend on x layer: B(x,y) = b_x * (-1)^y
        where a_y, b_x in {-1, 1} are random 1D arrays of size L.
        """
        # Random initial configuration for the chains
        a = np.random.choice([-1, 1], size=self.L)
        b = np.random.choice([-1, 1], size=self.L)
        
        A_spins = np.zeros((self.L, self.L), dtype=np.int8)
        B_spins = np.zeros((self.L, self.L), dtype=np.int8)
        
        for x in range(self.L):
            for y in range(self.L):
                A_spins[x, y] = a[y] * (1 if x % 2 == 0 else -1)
                B_spins[x, y] = b[x] * (1 if y % 2 == 0 else -1)
                
        return A_spins, B_spins

    def calc_energy(self, A_spins, B_spins):
        """
        Calculate total energy and energy per cell.
        """
        L = self.L
        # A-A bonds (horizontal x-direction)
        E_AA = np.sum(A_spins * np.roll(A_spins, -1, axis=0))
        # B-B bonds (vertical y-direction)
        E_BB = np.sum(B_spins * np.roll(B_spins, -1, axis=1))
        
        # A-B bonds (each A(x,y) connected to B(x,y), B(x-1,y), B(x,y-1), B(x-1,y-1))
        # Equivalently: sum over all x,y of A[x,y] * (B[x,y] + B[x-1,y] + B[x,y-1] + B[x-1,y-1])
        B_sum = B_spins + np.roll(B_spins, 1, axis=0) + np.roll(B_spins, 1, axis=1) + np.roll(np.roll(B_spins, 1, axis=0), 1, axis=1)
        E_AB = np.sum(A_spins * B_sum)
        
        E_total = E_AA + E_BB + E_AB
        return E_total, E_total / (L * L)

    def draw_configuration(self, A_spins, B_spins, filename):
        """
        Visualize the A and B spins.
        """
        fig, ax = plt.subplots(figsize=(6, 6))
        
        L = self.L
        # Plot A spins at integer coordinates
        XA, YA = np.meshgrid(np.arange(L), np.arange(L), indexing='ij')
        ax.scatter(XA.flatten(), YA.flatten(), c=A_spins.flatten(), cmap='bwr', s=100, marker='o', label='Sublattice A', edgecolors='k')
        
        # Plot B spins at half-integer coordinates
        XB, YB = np.meshgrid(np.arange(L) + 0.5, np.arange(L) + 0.5, indexing='ij')
        ax.scatter(XB.flatten(), YB.flatten(), c=B_spins.flatten(), cmap='bwr', s=100, marker='s', label='Sublattice B', edgecolors='k')
        
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1))
        plt.title(f'Ground State Configuration (L={L})')
        plt.tight_layout()
        os.makedirs('asset', exist_ok=True)
        plt.savefig(f'asset/{filename}', dpi=150)
        plt.close()

    def calc_correlation(self, num_samples=1000):
        """
        Sample multi-ground-states and calculate the correlation function C^{\mu \nu}(r).
        """
        L = self.L
        C_AA = np.zeros((L, L))
        C_BB = np.zeros((L, L))
        C_AB = np.zeros((L, L))
        
        for _ in range(num_samples):
            A, B = self.generate_ground_state()
            
            # Use FFT for fast correlation computation (convolution theorem)
            A_f = np.fft.fft2(A)
            B_f = np.fft.fft2(B)
            
            C_AA += np.real(np.fft.ifft2(A_f * np.conj(A_f))) / (L * L)
            C_BB += np.real(np.fft.ifft2(B_f * np.conj(B_f))) / (L * L)
            # For AB, r is the shift of B relative to A
            # B is physically at (x+0.5, y+0.5) but conceptually on the grid it's index offset
            C_AB += np.real(np.fft.ifft2(A_f * np.conj(B_f))) / (L * L)
            
        C_AA /= num_samples
        C_BB /= num_samples
        C_AB /= num_samples
        
        return C_AA, C_BB, C_AB

    def draw_correlation(self, corr, title, filename):
        fig, ax = plt.subplots(figsize=(6, 5))
        # Shift zero-frequency component to center
        # Since r takes integer values, let's just plot L x L
        # To show patterns well, we plot from -L/2 to L/2
        L = self.L
        C_shifted = np.fft.fftshift(corr)
        extent = [-L/2, L/2 - 1, -L/2, L/2 - 1]
        im = ax.imshow(C_shifted.T, origin='lower', cmap='RdBu', vmin=-1, vmax=1, extent=extent)
        plt.colorbar(im)
        ax.set_title(title)
        ax.set_xlabel('rx')
        ax.set_ylabel('ry')
        plt.tight_layout()
        plt.savefig(f'asset/{filename}', dpi=150)
        plt.close()

if __name__ == "__main__":
    np.random.seed(42)
    model = FrustratedIsing2D(L=16)
    
    print(f"--- 验证基态规则与能量 (L={model.L}) ---")
    E_theory = -2.0
    print(f"理论基态能量预测值：{E_theory} J (每元胞)")
    
    # 抽取 3 个样本并验证
    for i in range(3):
        A, B = model.generate_ground_state()
        E_tot, E_cell = model.calc_energy(A, B)
        print(f"样本 {i+1}: 观测总能量 = {E_tot:.1f}, 每元胞能量 = {E_cell:.3f} (与理论比对: {'正确' if np.isclose(E_cell, E_theory) else '异常'})")
        model.draw_configuration(A, B, f'sample_{i+1}.png')
        
    print("\n--- 计算并生成关联函数热力图 ---")
    print("开始进行蒙特卡洛随机基态采样，采样数量: 5000 ...")
    C_AA, C_BB, C_AB = model.calc_correlation(num_samples=5000)
    
    model.draw_correlation(C_AA, 'Correlation $C^{AA}(r)$', 'corr_AA.png')
    model.draw_correlation(C_BB, 'Correlation $C^{BB}(r)$', 'corr_BB.png')
    
    print("关联函数已处理完毕并保存在 asset/ 目录下。")
    
    print("\n观察：")
    print("C_AA 对角线方向 (rx != 0, ry != 0):", C_AA[1, 1], C_AA[2, 2])
    print("C_BB 对角线方向 (rx != 0, ry != 0):", C_BB[1, 1], C_BB[2, 2])
    print("A 和 B 之间的交叉关联 C_AB:", np.max(np.abs(C_AB)))
