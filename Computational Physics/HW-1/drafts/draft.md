
`basic_component.py`：

```python
import numpy as np

class GridPlaygroundContainer():
	def __init__(self, L):
		self.L = L
		# self.t_max=t_max
		#self.init_playground = np.array(np.zeros(L,L))
		#return self.init_playground
		
	def add_cube(self, playground: array):
		'''
		输入：正方形数组
		操作：随机挑一个位置+1
		'''
		pass
		
	def check_playground(self, playground: array):
		'''
		输入：正方形数组
		操作：并行检查一次是否有大于等于4的位置，返回一个布尔值。大于等于4的位置则为True
		注意仅检查一次！
		'''
		pass
		
	def minus_four(self, 
				playground: array, 
				check_array: array):
		'''
		输入：布尔值的正方形数组; playground正方形数组
		操作：为True的位置，playground的数组的值-4
		返回：playground
		'''
		pass
		
	def cube_transfer(self,
					playground: array, 
					check_array: array):
		'''
		输入：布尔值的正方形数组; playground正方形数组
		操作：每个True的位置，其上下左右位置playground数值+1
		需注意位于边界的点，边界外的不用处理
		返回
		'''
		pass
	
	def sum_cube(self,  playground:array):
		'''
		计算所有的cube数量
		'''
		return np.sum(playground)
	
	def one_step(self, playground:array):
		'''
		注意此处仅为示意，语法可能有错
		
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
		init_check = sellf.check_playground(playground)
		if init_check == True:
			raise()
			
		playground = self.add_cube(playground)
		
		score = 0
		while True:
			check_array = self.check_playground(playground)
			score += np.sum(check_array)
			# 如果没有True的项，则结束
			if check_array == False:
				break
			
			playground = self.minus_four(playground, check_array)
			playground = self.cube_transfer(playground, check_array)
			
			
		return playground, score, sum_cube(playground)/self.L**2
		
```



`question_1_trail.py`:

```python
# import ...

if __name__ == "__main__":
	L =32
	container = GridPlaygroundContainer(L)
	t_steps = 10000
	
	# 初始化
	playground = np.array(np.zeros(L,L))
	# scores = np.zeros(t_steps)
	cube_densities = np.zeros(t_steps)
	
	for t_step in range(t_steps):
		playground, _, cube_density = container.one_step(playgorund)
		# scores[t_step] = score
		cube_densities[t_step+1] = cube_density
```

`question_2_trail.py`:

```python
# import ...

def get_distribution(scores: array):
	'''
	输入：记录分数的数组
	输出：F为每个分数对应的频率数组
	注意计算的高效性
	'''
	return F

if __name__ == "__main__":
	L =64
	container = GridPlaygroundContainer(L)
	t_steps = 10000
	
	# 初始化
	playground = np.array(np.zeros(L,L))
	scores = np.zeros(t_steps)
	# cube_densities = np.zeros(t_steps)
	
	for t_step in range(t_steps):
		playground, score, _ = container.one_step(playgorund)
		scores[t_step] = score
		# cube_densities[t_step+1] = cube_density
		
	F = get_distribution(scores)
	# 画图
```

`container_parallel.py`:

```python
# import ...

class ContainerParallel():
	'''
	并行地实例化不同尺寸的容器
	并计算其分布，并提供绘图
	'''
	def __init__(self, Ls: list, t_step_max: int = 5000):
		'''
		Ls为不同的尺寸的list
		t_step_max为时间步
		'''
		pass
		
	def get_distribution_par(self, scores: array):
		'''
		功能和前面的一样,但是是并行的
		输出也为array
		'''
		pass
		
	def plot_distribution(self, F: array):
		'''
		输入为欸上面得到的多维数组, 画在一个表中
		'''
		pass
	
```
