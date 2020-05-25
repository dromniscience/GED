# -*- coding:utf-8 -*-

"""
lasp.py: LSAP Solver

    Currently this file contains a Munkres' Hungarian implementation.

Written by Ding Rui
Latest Version: 2020/5/25
"""
import numpy as np
import time


class Munkres:
    def __init__(self):
        pass

    def __call__(self, matrix: np.ndarray, maximize=False) -> (np.ndarray, np.ndarray):
        # 要求矩阵是2维方阵, 并且每个元素都是整数
        self.CheckSafety(matrix)
        # 构造实例变量
        self.len = int(matrix.shape[0])
        self.matrix = (-matrix if maximize else matrix.copy())  # 总是修改原始数据的一个副本
        self.colcover = np.full(self.len, False)
        self.rowcover = np.full(self.len, False)
        self.step = 0
        self.zprime = set()
        self.zstar = set()

        temp = None  # 用于记忆 Step 1 可能给出的 0' 位置
        while True:
            if self.step == 0:
                self.Step0()
            elif self.step == 1:
                temp = self.Step1()
            elif self.step == 2:
                self.Step2(temp)
            elif self.step == 3:
                self.Step3()
            else:
                break

        # 由 0* 构造指派
        row = np.arange(self.len)
        col = np.array([j for _, j in sorted(list(self.zstar))])
        return row, col

    # 检查输入是否合法
    def CheckSafety(self, matrix: np.ndarray):
        if len(matrix.shape) != 2:
            raise ValueError("expected a matrix (2-d array), got a %r array"
                             % (matrix.shape,))
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("expected a square matrix (2-d array), got a %r array"
                             % (matrix.shape,))
        if not np.issubdtype(matrix.dtype, np.integer):
            raise ValueError("expected a matrix containing integer entries, got %s"
                             % (matrix.dtype,))

    # 找到初始的一组极大的独立零元素并标 0*
    def InitZeros(self):
        # 所有零元素的集合
        zeros = set(map(tuple, np.argwhere(self.matrix == 0).tolist()))

        while zeros:
            zero = zeros.pop()
            self.zstar.add(zero)
            for i in range(self.len):
                if (i, zero[1]) in zeros:
                    zeros.remove((i, zero[1]))
                if (zero[0], i) in zeros:
                    zeros.remove((zero[0], i))

    # 寻找一个尚未覆盖的 0, 没有则返回 ()
    def UncoveredZero(self):
        temp = self.matrix[~self.rowcover, :]
        temp = temp[:, ~self.colcover]
        temp = np.argwhere(temp == 0)
        if int(temp.shape[0]) == 0:  # 全都覆盖了
            return ()
        else:
            zero = tuple(temp[0].tolist())
            # 确定在原始图的坐标
            temp, x, y = -1, 0, 0
            for row in range(self.len):
                if self.rowcover[row] == False:
                    temp += 1
                if temp == zero[0]:
                    x = row
                    break
            temp = -1
            for col in range(self.len):
                if self.colcover[col] == False:
                    temp += 1
                if temp == zero[1]:
                    y = col
                    break
            return (x, y)

    # 找出 0' 所在列上的那个 0* (可能存在)
    def FindColStar(self, zero) -> tuple:
        for i in range(self.len):
            if (i, zero[1]) in self.zstar:
                return (i, zero[1])
        return ()

    # 找出 0* 所在行上的那个 0' (必然存在)
    def FindRowPrime(self, zero) -> tuple:
        for i in range(self.len):
            if (zero[0], i) in self.zprime:
                return (zero[0], i)

    # 初始化
    def Step0(self):
        # 每列都减去相应最小值, 然后每行类似
        self.matrix -= self.matrix.min(axis=0, keepdims=True)  # type: np.ndarray
        self.matrix -= self.matrix.min(axis=1, keepdims=True)  # type: np.ndarray
        # 找到一组极大独立零元素, 标记为 0*
        self.InitZeros()
        # 列覆盖 0*
        for zero in self.zstar:
            self.colcover[zero[1]] = True
        # 结束此步: 如果所有列均覆盖, 自然已找到解; 否则去 Step 1
        self.step = 1 if self.colcover.sum() < self.len else 4
        return

    # Step 1
    def Step1(self):
        # 找到一个未覆盖的 0
        zero = self.UncoveredZero()
        if not zero:  # 未覆盖的全是正整数, 可以等价变换代价矩阵 Step 3
            self.step = 3
            return
        # 标记为 0'
        self.zprime.add(zero)
        for i in range(self.len):
            if (zero[0], i) in self.zstar:  # 同行有 0*, 将对应列覆盖改成行覆盖
                self.colcover[i] = False
                self.rowcover[zero[0]] = True
                self.step = 1
                return
        # 0' 所在行无 0*, 说明这时有增广路径
        self.step = 2
        return zero

    # 增广独立零元素的数量
    def Step2(self, zero: tuple):
        # 找增广路径
        path = {zero}
        zero = self.FindColStar(zero)
        while zero:
            path.add(zero)
            zero = self.FindRowPrime(zero)
            path.add(zero)
            zero = self.FindColStar(zero)
        # 更新标记
        # 将增广路上的 0' 改成 0*, 0* 标记去掉; 这时 0* 数量恰好比之前多一个
        self.zstar = self.zstar ^ path
        # 清空 0' 标记
        self.zprime = set()
        # 重新覆盖
        self.rowcover = np.full(self.len, False)
        self.colcover = np.full(self.len, False)
        for zero in self.zstar:
            self.colcover[zero[1]] = True
        # 结束此步: 如果所有列均覆盖, 自然已找到解; 否则去 Step 1
        self.step = 1 if self.colcover.sum() < self.len else 4
        return

    # 等价变换代价矩阵
    def Step3(self):
        # 寻找未覆盖的元素中的最小值 (一定是正整数)
        temp = self.matrix[~self.rowcover, :]  # type: np.ndarray
        temp = temp[:, ~self.colcover]
        decrease = temp.min()

        for i in range(self.len):
            if self.rowcover[i] == False:
                self.matrix[i, :] -= decrease
        for i in range(self.len):
            if self.colcover[i] == True:
                self.matrix[:, i] += decrease

        # 其他参数无需更新
        self.step = 1
        return


if __name__ == '__main__':
    import random

    random.seed(1)

    from scipy.optimize import linear_sum_assignment

    myLSAP = Munkres()

    data = []
    for _ in range(100):
        scale = 100
        data.append(np.random.randint(-100, 100, (scale, scale)))

    ansMy = []
    ansScipy = []

    st = time.time()
    for mat in data:
        row, col = myLSAP(mat)
        ansMy.append(mat[row, col].sum())
    en = time.time()
    mytime = en - st

    st = time.time()
    for mat in data:
        row, col = linear_sum_assignment(mat)
        ansScipy.append(mat[row, col].sum())
    en = time.time()
    Scipytime = en - st

    print(ansMy == ansScipy)
    print(mytime / Scipytime)
