# -*- coding: utf-8 -*-
"""
TSSP方法求解一维含时Gross-Pitaevskii方程
Time-Splitting Spectral Method for 1D Gross-Pitaevskii Equation

作者: kyksj-1
物理背景: 薛定谔视角下的量子波包动力学

方程:
    i * dψ/dt = -1/2 * d²ψ/dx² + V(x)ψ + η(ψ)ψ
    其中 V(x) = x²/2 (谐振子势), η(ψ) = |ψ|²/2 (非线性项)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

# 设置中文字体
rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
rcParams['axes.unicode_minus'] = False


class TSSPGSolver:
    """
    TSSP方法求解GP方程的核心求解器
    
    算法原理:
        将哈密顿量H分解为动能项T和势能项(含非线性项)V:
        H = T + V, 其中 T = -1/2 * d²/dx², V = V(x) + |ψ|²/2
        
        时间演化算子 exp(-iHΔt) 通过Strang分裂近似:
        exp(-iHΔt) ≈ exp(-iVΔt/2) * exp(-iTΔt) * exp(-iVΔt/2)
        
        动能项在动量空间处理(谱方法)，势能项在坐标空间处理
    """
    
    def __init__(self, L=20.0, N=512, dt=0.01, t_max=20.0):
        """
        初始化求解器
        
        参数:
            L: 空间范围 [-L, L]
            N: 空间网格点数
            dt: 时间步长
            t_max: 最大演化时间
        """
        self.L = L
        self.N = N
        self.dt = dt
        self.t_max = t_max
        
        # 空间网格
        self.x = np.linspace(-L, L, N, endpoint=False)
        self.dx = 2 * L / N
        
        # 动量网格 (用于谱方法)
        self.k = np.fft.fftfreq(N, d=self.dx) * 2 * np.pi
        
        # 势能项 V(x) = x²/2
        self.V = 0.5 * self.x**2
        
        # 时间步数
        self.num_steps = int(t_max / dt)
        self.t = np.linspace(0, t_max, self.num_steps + 1)
        
        # 存储结果
        self.rho_history = []  # 密度演化历史
        self.width_history = []  # 波包宽度演化历史
        
    def initialize_wavefunction(self):
        """
        初始化波函数
        |ψ(x,0)| = (1/sqrt(2π)) * exp(-x²/2)
        选择相位为零的实函数
        """
        psi0 = np.exp(-self.x**2 / 2) / np.sqrt(2 * np.pi)
        return psi0.astype(np.complex128)
    
    def kinetic_step(self, psi, dt):
        """
        动能项演化 (在动量空间处理)
        exp(-i * T * dt) * ψ = IFFT[exp(i * k² * dt/2) * FFT[ψ]]
        
        物理意义: 自由粒子的演化，波包扩散
        """
        psi_k = np.fft.fft(psi)
        psi_k *= np.exp(-0.5j * self.k**2 * dt)
        return np.fft.ifft(psi_k)
    
    def potential_step(self, psi, dt):
        """
        势能项演化 (含非线性项，在坐标空间处理)
        exp(-i * V_eff * dt) * ψ
        其中 V_eff = V(x) + |ψ|²/2
        
        物理意义: 局域相位演化，无空间扩散
        """
        V_eff = self.V + 0.5 * np.abs(psi)**2
        psi *= np.exp(-1j * V_eff * dt)
        return psi
    
    def tssp_step(self, psi):
        """
        单步TSSP演化 (Strang分裂)
        ψ(t+dt) = exp(-iV*dt/2) * exp(-iT*dt) * exp(-iV*dt/2) * ψ(t)
        
        复杂度分析:
        - 每步需要2次FFT和2次IFFT: O(N log N)
        - 总时间复杂度: O(M * N log N), M为时间步数
        """
        psi = self.potential_step(psi, self.dt / 2)
        psi = self.kinetic_step(psi, self.dt)
        psi = self.potential_step(psi, self.dt / 2)
        return psi
    
    def compute_density(self, psi):
        """计算密度函数 ρ(x) = |ψ(x)|²"""
        return np.abs(psi)**2
    
    def compute_width(self, psi):
        """
        计算波包宽度 w(t) = <x²> = ∫ x²|ψ|² dx
        物理意义: 波包的空间展宽程度
        """
        rho = self.compute_density(psi)
        return np.sum(self.x**2 * rho) * self.dx
    
    def evolve(self, save_interval=10):
        """
        时间演化主循环
        
        参数:
            save_interval: 每隔多少步保存一次结果
        """
        psi = self.initialize_wavefunction()
        
        print("=" * 60)
        print("TSSP方法求解GP方程 - 时间演化开始")
        print("=" * 60)
        print(f"空间范围: [-{self.L}, {self.L}]")
        print(f"网格点数: {self.N}")
        print(f"时间步长: {self.dt}")
        print(f"演化时间: [0, {self.t_max}]")
        print(f"总时间步数: {self.num_steps}")
        print("=" * 60)
        
        # 保存初始状态
        self.rho_history.append(self.compute_density(psi).copy())
        self.width_history.append(self.compute_width(psi))
        
        # 时间演化
        for step in range(1, self.num_steps + 1):
            psi = self.tssp_step(psi)
            
            # 归一化校正 (保持概率守恒)
            norm = np.sqrt(np.sum(np.abs(psi)**2) * self.dx)
            psi /= norm
            
            if step % save_interval == 0:
                self.rho_history.append(self.compute_density(psi).copy())
                self.width_history.append(self.compute_width(psi))
                
                if step % (self.num_steps // 10) == 0:
                    t_current = step * self.dt
                    width = self.width_history[-1]
                    print(f"t = {t_current:6.2f}, 波包宽度 = {width:.6f}")
        
        print("=" * 60)
        print("演化完成!")
        print("=" * 60)
        
        return psi
    
    def plot_density_heatmap(self, save_path='../asset/density_heatmap.png'):
        """
        绘制密度演化热力图
        """
        rho_array = np.array(self.rho_history).T
        t_saved = np.linspace(0, self.t_max, len(self.rho_history))
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 选择合适的显示范围
        x_idx = np.abs(self.x) <= 8
        im = ax.imshow(rho_array[x_idx, :], 
                       extent=[0, self.t_max, -8, 8],
                       aspect='auto', origin='lower',
                       cmap='hot', vmin=0)
        
        cbar = plt.colorbar(im, ax=ax, label='密度 ρ(x,t)')
        ax.set_xlabel('时间 t', fontsize=14)
        ax.set_ylabel('空间坐标 x', fontsize=14)
        ax.set_title('密度函数 ρ(x,t) 演化热力图\n(GP方程, TSSP方法)', fontsize=16)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"热力图已保存至: {save_path}")
        
    def plot_width_evolution(self, save_path='../asset/width_evolution.png'):
        """
        绘制波包宽度演化
        """
        t_saved = np.linspace(0, self.t_max, len(self.width_history))
        width_array = np.array(self.width_history)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(t_saved, width_array, 'b-', linewidth=2, label='数值解')
        ax.axhline(y=width_array[0], color='r', linestyle='--', 
                   label=f'初始宽度 = {width_array[0]:.4f}')
        
        ax.set_xlabel('时间 t', fontsize=14)
        ax.set_ylabel('波包宽度 w(t) = ⟨x²⟩', fontsize=14)
        ax.set_title('波包宽度随时间演化\n(GP方程, TSSP方法)', fontsize=16)
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"宽度演化图已保存至: {save_path}")
        
    def analyze_results(self):
        """
        分析数值结果
        """
        width_array = np.array(self.width_history)
        
        print("\n" + "=" * 60)
        print("数值结果分析")
        print("=" * 60)
        print(f"初始波包宽度 w(0) = {width_array[0]:.6f}")
        print(f"最终波包宽度 w({self.t_max}) = {width_array[-1]:.6f}")
        print(f"平均波包宽度 = {np.mean(width_array):.6f}")
        print(f"宽度标准差 = {np.std(width_array):.6f}")
        print(f"宽度最大值 = {np.max(width_array):.6f}")
        print(f"宽度最小值 = {np.min(width_array):.6f}")
        
        # 检查周期性
        if len(width_array) > 100:
            from scipy.signal import find_peaks
            peaks, _ = find_peaks(width_array, height=np.mean(width_array))
            if len(peaks) > 1:
                period = np.mean(np.diff(peaks)) * (self.t_max / len(width_array))
                print(f"检测到振荡周期 ≈ {period:.4f}")
        
        print("=" * 60)


def main():
    """主函数"""
    solver = TSSPGSolver(L=20.0, N=512, dt=0.01, t_max=20.0)
    solver.evolve(save_interval=10)
    
    solver.plot_density_heatmap()
    solver.plot_width_evolution()
    solver.analyze_results()
    
    return solver


if __name__ == '__main__':
    solver = main()
