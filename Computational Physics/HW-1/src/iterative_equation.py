'''
迭代映射及其导数
Q1: x_{n+1} = 1 - mu * x_n^2
Q2: x_{n+1} = cos(x_n) - mu * x_n^2
'''
import numpy as np


def linear_iter(x, mu=0.5):
    return 1 - mu * x**2

def linear_deriv(x, mu=0.5):
    return -2 * mu * x

def cosine_iter(x, mu=0.5):
    return np.cos(x) - mu * x**2

def cosine_deriv(x, mu=0.5):
    return -np.sin(x) - 2 * mu * x
