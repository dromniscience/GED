# -*- coding:utf-8 -*-

"""
IPFP1.py: QAP using IPFP approximation strategy
    
    Initialize with the algorithm identical to the one in BGM1.py

Written by Ding Rui
Last Version: 2020/5/24

"""

import numpy as np
from common import ArgParse, Graph, CostMatrix, PrintGED, Add1Star, SolveLSAP, DeltaMatrix

dot_sub, dot_ins, dot_del, edge_sub, edge_ins, edge_del, root_path, void, inf,\
    g1, g2 = ArgParse()

# 1-Star 的 LSAP 初始解
def LSAP1Star(graph1: Graph, graph2: Graph):
    return SolveLSAP(CostMatrix(graph1, graph2, Add1Star))

# IPFP 的算法实现
def IPFP(X: np.ndarray, cost: np.ndarray, delta: np.ndarray, iter=40):
    Linear = (cost*X).sum()
    Sum = 0.5 * (X.reshape(-1) @ delta @ X.reshape(-1)) + Linear

    while iter > 0:
        iter -= 1
        # 在 X 的 Jacobbi 矩阵  Jacobbi = (tot, tot)
        Jacobbi = (X.reshape(-1) @ delta).reshape(X.shape) + cost
        # 梯度最大位置的反向  B_ = (tot, tot)
        _, _, B_ = SolveLSAP(Jacobbi)
        # 新的位置的目标值
        Linear_ = (cost * B_).sum()
        Sum_ = 0.5 * (B_.reshape(-1) @ delta @ B_.reshape(-1)) + Linear_
        # 求极值点位置
        alpha = (Jacobbi * B_).sum() - 2 * Sum + Linear
        beta = Sum_ + Sum - (Jacobbi * B_).sum() - Linear
        beta = (beta if abs(beta) >= 1e-8 else -1) # 防止除零, 此操作不改变分支走向
        extrema = - alpha / (2 * beta)

        # 更新参数
        temp = X   # 临时记忆, 观察此次有无真正更新
        if beta <= 0 or extrema >= 1:
            X, Linear, Sum = B_, Linear_, Sum_
        else:
            X = X + extrema * (B_ - X)
            Linear = (cost * X).sum()
            Sum = Sum - alpha ** 2 / (4 * beta)
        # 已找到局部极值
        if abs(X - temp).sum() <= 1e-5:
            break
    # 寻找最接近的 map
    # print(iter)
    _, col, _ = SolveLSAP(-X)
    return col

g1 = Graph(root_path + '/'+ g1)
g2 = Graph(root_path + '/'+ g2)

cost, delta = DeltaMatrix(g1, g2)
_, _, init = LSAP1Star(g1, g2)

col = IPFP(init, cost, delta)
answer = [(int(col[i]) if col[i] < g2.dots else void) for i in range(g1.dots)]
PrintGED(g1, g2, tuple(answer))
