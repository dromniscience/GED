# -*- coding:utf-8 -*-

"""
IPFP0.py: QAP using IPFP approximation strategy

    Initialize with the algorithm identical to the one in BGM1.py

Written by Ding Rui
Version 1.0: 2020/5/24
Version 1.1: 2020/5/27  Compulsory exit when running time reaches 30s
Version 1.2: 2020/5/29  Reduced matrix applied, speed up the process
"""

import numpy as np
import time
from reduced import inf, ArgParse, Graph, CostMatrix, PrintGED, Add1Star, SolveLSAP, ReducedDeltaMatrix, Padding, Unpadding

dot_sub, dot_ins, dot_del, edge_sub, edge_ins, edge_del, root_path, void, inf,\
    g1, g2 = ArgParse()

# 计时开始
start = time.time()

# 1-Star 的 LSAP 初始解
def ReducedLSAP1Star(graph1: Graph, graph2: Graph):
    return Unpadding(SolveLSAP(CostMatrix(graph1, graph2, Add1Star))[-1], graph1.dots, graph2.dots)

# IPFP 的算法实现
def ReducedIPFP(X: np.ndarray, cost: np.ndarray, delta: np.ndarray, iter=40):
    global start, inf

    Linear = (cost*X).sum()
    Sum = 0.5 * (X.reshape(-1) @ delta @ X.reshape(-1)) + Linear

    while iter > 0:
        iter -= 1

        # 强制退出
        if time.time() - start >= 29:
            break

        # 在 X 的 Jacobbi 矩阵  Jacobbi = (tot, tot)
        Jacobbi = (X.reshape(-1) @ delta).reshape(X.shape) + cost
        # 梯度最大位置的反向
        B_ = Unpadding(SolveLSAP(Padding(Jacobbi, inf))[-1], X.shape[0]-1, X.shape[1]-1)
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
    _, col, _ = SolveLSAP(Padding(-X, 0))
    return col

g1 = Graph(root_path + '/'+ g1)
g2 = Graph(root_path + '/'+ g2)

cost, delta = ReducedDeltaMatrix(g1, g2)
col = ReducedIPFP(ReducedLSAP1Star(g1, g2), cost, delta)
answer = [(int(col[i]) if col[i] < g2.dots else void) for i in range(g1.dots)]
PrintGED(g1, g2, tuple(answer))
