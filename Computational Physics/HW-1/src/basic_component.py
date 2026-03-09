# basic_component.py
import numpy as np

class GridPlaygroundContainer():
    def __init__(self, L):
        self.L = L

    def add_cube(self, playground: np.ndarray):
        '''
        输入：正方形数组
        操作：随机挑一个位置+1
        '''
        i, j = np.random.randint(0, self.L, size=2)
        playground[i, j] += 1
        return playground

    def check_playground(self, playground: np.ndarray):
        '''
        输入：正方形数组
        操作：并行检查一次是否有大于等于4的位置，返回一个布尔值。大于等于4的位置则为True
        注意仅检查一次！
        '''
        return playground >= 4

    def minus_four(
        self,
        playground: np.ndarray,
        check_array: np.ndarray,
    ):
        '''
        输入：布尔值的正方形数组; playground正方形数组
        操作：为True的位置，playground的数组的值-4
        返回：playground
        '''
        playground[check_array] -= 4
        return playground

    def cube_transfer(
        self,
        playground: np.ndarray,
        check_array: np.ndarray,
    ):
        '''
        输入：布尔值的正方形数组; playground正方形数组
        操作：每个True的位置，其上下左右位置playground数值+1
        需注意位于边界的点，边界外的不用处理
        返回：playground
        具体实现：
        利用错位切片进行四周溢出。
        越界部分被切片机制自然丢弃，完美模拟开放边界条件的耗散。
        '''
        topples = check_array.astype(int)

        playground[:-1, :] += topples[1:, :]
        playground[1:, :] += topples[:-1, :]
        playground[:, :-1] += topples[:, 1:]
        playground[:, 1:] += topples[:, :-1]

        return playground

    def sum_cube(self, playground: np.ndarray):
        '''
        计算所有的cube数量
        '''
        return np.sum(playground)

    def one_step(self, playground: np.ndarray):
        '''
        完成一步的迭代：
        输入一个playground
        先检查是否是符合要求的（完成了上一步的）
        然后放上去一个cube
        然后开始迭代检查：
        先把满4个的给清了；
        然后加到它们前后左右.
        这就完成了一次迭代，等它再次检查，是否会break

        返回：这一步完成的playground、得分和方块密度
        '''
        init_check = self.check_playground(playground)
        if np.any(init_check):
            raise ValueError('初始状态不符合要求')

        playground = self.add_cube(playground)

        score = 0
        while True:
            check_array = self.check_playground(playground)
            if not np.any(check_array):
                break
            score += int(np.sum(check_array))
            playground = self.minus_four(playground, check_array)
            playground = self.cube_transfer(playground, check_array)

        return playground, score, self.sum_cube(playground) / (self.L**2)
