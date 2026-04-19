# -*- coding: utf-8 -*-
"""
堆上的最短路径问题求解
使用动态规划算法找最短路径，分析统计规律

Author: kyksj-1
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List
import os

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False


class HeapStructure:
    """
    N层堆结构
    
    第n层有n+1个节点，每个节点值为[0,1)均匀分布随机数
    节点v_i^n指向下一层的两个子节点v_i^{n+1}和v_{i+1}^{n+1}
    """
    
    def __init__(self, N: int, seed: int = None):
        """
        初始化N层堆
        
        Args:
            N: 堆的层数（从0到N）
            seed: 随机种子（用于可重复性）
        """
        self.N = N
        if seed is not None:
            np.random.seed(seed)
        
        self.values = []
        for n in range(N + 1):
            layer_values = np.random.uniform(0, 1, n + 1)
            self.values.append(layer_values)
    
    def find_shortest_path_dp(self) -> Tuple[List[Tuple[int, int]], float, int]:
        """
        使用动态规划找最短路径
        
        从顶部到底部，dp[n][i]表示到达第n层第i个节点的最小路径长度
        dp[n+1][i] = min(dp[n][i-1], dp[n][i]) + values[n+1][i]
        dp[n+1][i+1] = min(dp[n][i-1], dp[n][i]) + values[n+1][i+1]
        
        Returns:
            path: 最短路径，每个元素为(层号, 节点索引)
            length: 最短路径长度
            x_star: 终点横坐标
        """
        dp = []
        parent = []
        
        dp.append([self.values[0][0]])
        parent.append([(-1, -1)])
        
        for n in range(1, self.N + 1):
            dp_layer = []
            parent_layer = []
            
            for i in range(n + 1):
                candidates = []
                
                if i > 0:
                    candidates.append((dp[n-1][i-1], (n-1, i-1)))
                
                if i < n:
                    candidates.append((dp[n-1][i], (n-1, i)))
                
                min_val, min_parent = min(candidates, key=lambda x: x[0])
                
                dp_layer.append(min_val + self.values[n][i])
                parent_layer.append(min_parent)
            
            dp.append(dp_layer)
            parent.append(parent_layer)
        
        min_length = min(dp[self.N])
        final_idx = dp[self.N].index(min_length)
        
        path = [(self.N, final_idx)]
        n, i = self.N, final_idx
        while n > 0:
            i = parent[n][i][1]
            n -= 1
            path.append((n, i))
        
        path.reverse()
        
        return path, min_length, final_idx
    
    def get_path_values(self, path: List[Tuple[int, int]]) -> List[float]:
        """获取路径上所有节点的值"""
        return [self.values[n][i] for n, i in path]


def analyze_w_vs_N(N_range: range, num_samples: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
    """
    分析w(N)随N的变化规律
    
    w(N) = sqrt(<[x*(N)]^2> - <x*(N)>^2)
    
    Args:
        N_range: N的取值范围
        num_samples: 每个N的样本数
    
    Returns:
        N_array: N值数组
        w_array: w(N)值数组
    """
    N_list = list(N_range)
    w_results = []
    
    for N in N_list:
        x_stars = []
        
        for sample in range(num_samples):
            heap = HeapStructure(N, seed=None)
            _, _, x_star = heap.find_shortest_path_dp()
            x_stars.append(x_star)
        
        x_stars = np.array(x_stars)
        
        mean_x = np.mean(x_stars)
        mean_x2 = np.mean(x_stars ** 2)
        w = np.sqrt(mean_x2 - mean_x ** 2)
        
        w_results.append(w)
        print(f"N = {N:3d}, w(N) = {w:.4f}, <x*> = {mean_x:.4f}, 样本数 = {num_samples}")
    
    return np.array(N_list), np.array(w_results)


def plot_w_vs_N(N_array: np.ndarray, w_array: np.ndarray, save_path: str = None):
    """绘制w(N)随N的变化"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].plot(N_array, w_array, 'o-', markersize=6, linewidth=2, color='#2E86AB')
    axes[0].set_xlabel('N (堆层数)', fontsize=12)
    axes[0].set_ylabel('w(N)', fontsize=12)
    axes[0].set_title('w(N)随N的变化规律', fontsize=14)
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(N_array, w_array, 'o-', markersize=6, linewidth=2, color='#2E86AB', label='数值结果')
    
    if len(N_array) > 2:
        coeffs = np.polyfit(N_array, w_array, 1)
        fit_line = np.polyval(coeffs, N_array)
        axes[1].plot(N_array, fit_line, '--', linewidth=2, color='#E94F37', 
                     label=f'线性拟合: w = {coeffs[0]:.3f}N + {coeffs[1]:.3f}')
        axes[1].legend(fontsize=10)
    
    axes[1].set_xlabel('N (堆层数)', fontsize=12)
    axes[1].set_ylabel('w(N)', fontsize=12)
    axes[1].set_title('w(N)与线性拟合对比', fontsize=14)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"w(N)变化图已保存至: {save_path}")
    
    plt.close()


def plot_single_heap_example(N: int = 5, save_path: str = None):
    """绘制单个堆结构示例及最短路径"""
    heap = HeapStructure(N, seed=42)
    path, length, x_star = heap.find_shortest_path_dp()
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for n in range(N + 1):
        for i in range(n + 1):
            x = i - n / 2
            y = N - n
            ax.scatter(x, y, s=500, c='#2E86AB', zorder=3)
            ax.text(x, y, f'{heap.values[n][i]:.2f}', ha='center', va='center', 
                   fontsize=8, color='white', fontweight='bold', zorder=4)
    
    for n in range(N):
        for i in range(n + 1):
            x_parent = i - n / 2
            y_parent = N - n
            y_child = N - n - 1
            
            x_child1 = i - (n + 1) / 2
            x_child2 = i + 1 - (n + 1) / 2
            
            ax.plot([x_parent, x_child1], [y_parent, y_child], 'k-', linewidth=1, alpha=0.3)
            ax.plot([x_parent, x_child2], [y_parent, y_child], 'k-', linewidth=1, alpha=0.3)
    
    if len(path) > 1:
        path_x = []
        path_y = []
        for n, i in path:
            path_x.append(i - n / 2)
            path_y.append(N - n)
        ax.plot(path_x, path_y, 'r-', linewidth=3, zorder=5, label='最短路径')
        ax.scatter(path_x, path_y, s=600, c='#E94F37', zorder=6)
    
    ax.set_xlim(-N/2 - 1, N/2 + 1)
    ax.set_ylim(-0.5, N + 0.5)
    ax.set_xlabel('横坐标 (i - n/2)', fontsize=12)
    ax.set_ylabel('层号 (从下往上)', fontsize=12)
    ax.set_title(f'N={N}层堆结构示例\n最短路径长度={length:.4f}, 终点x*={x_star}', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2)
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"堆结构示例图已保存至: {save_path}")
    
    plt.close()
    
    return heap, path, length, x_star


def main():
    """主函数"""
    print("=" * 60)
    print("堆上的最短路径问题求解")
    print("=" * 60)
    
    asset_dir = os.path.join(os.path.dirname(__file__), '..', 'asset')
    os.makedirs(asset_dir, exist_ok=True)
    
    print("\n" + "=" * 60)
    print("1. 单个堆结构示例")
    print("=" * 60)
    
    heap, path, length, x_star = plot_single_heap_example(
        N=6, 
        save_path=os.path.join(asset_dir, 'heap_example.png')
    )
    
    print(f"\n堆层数 N = 6")
    print(f"最短路径: {path}")
    print(f"最短路径长度: {length:.4f}")
    print(f"终点横坐标 x*: {x_star}")
    print(f"路径上节点值: {[f'{v:.4f}' for v in heap.get_path_values(path)]}")
    
    print("\n" + "=" * 60)
    print("2. w(N)随N的变化规律分析")
    print("=" * 60)
    print("对每个N生成1000个不同的堆，计算w(N)")
    print()
    
    N_array, w_array = analyze_w_vs_N(range(5, 31, 5), num_samples=1000)
    
    plot_w_vs_N(N_array, w_array, save_path=os.path.join(asset_dir, 'w_vs_N.png'))
    
    print("\n" + "=" * 60)
    print("数值结果汇总")
    print("=" * 60)
    print(f"N值范围: {N_array[0]} 到 {N_array[-1]}")
    print(f"w(N)最小值: {w_array.min():.4f} (N={N_array[w_array.argmin()]})")
    print(f"w(N)最大值: {w_array.max():.4f} (N={N_array[w_array.argmax()]})")
    
    if len(N_array) > 2:
        coeffs = np.polyfit(N_array, w_array, 1)
        print(f"\n线性拟合结果: w(N) = {coeffs[0]:.4f} * N + {coeffs[1]:.4f}")
        print(f"拟合斜率 (增长率): {coeffs[0]:.4f}")
    
    print("\n" + "=" * 60)
    print("算法复杂度分析")
    print("=" * 60)
    print("动态规划算法:")
    print("  - 总节点数: O(N^2)")
    print("  - 每个节点处理: O(1)")
    print("  - 总时间复杂度: O(N^2)")
    print("  - 空间复杂度: O(N^2) (存储节点值和DP表)")


if __name__ == '__main__':
    main()
