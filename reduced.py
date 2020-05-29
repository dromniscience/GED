# -*- coding:utf-8 -*-

'''
common.py: includes the definition of reusable codes, namely:

        ArgParse    funct   Collecting & Parsing command line arguments
        Graph       class   Reorganize data from input files
        ComputeGED  funct   Calculate GED given a dot-wise mapping
        CostMatrix  funct   Compute for cost matrix given a relevant function
        PrintGED    funct   Generate Formatted output
        Add1Star    funct   Compute 1-Star local cost
        SolveLSAP   funct   A LSAP solver
        DeltaMatrix funct   Formulate first-order and second-order cost matrix

    Warning: Please do NOT modify those values returned by ArgParse!
             Otherwise it may produce incorrect answer for GED.

Written by Ding Rui
Latest Version: 2020/5/24
'''

# from Common import * 时仅允许导入以下变量
# 但建议不要使用此格式导入本模块!

import argparse
import numpy as np
from scipy.optimize import linear_sum_assignment

root_path = r'./GED_data/preprocessed_C/'
void = None
inf = 0x7FFFFFFF  # 当然更大也行。这是一个潜在的风险，更好的做法是由输入推断设置上界。

dot_sub = None
dot_ins = None
dot_del = None
edge_sub = None
edge_ins = None
edge_del = None


# 差异集计算
def InplaceDifferDict(dict1: dict, dict2: dict):
    for key in dict2:
        dict1[key] = dict1.get(key, 0) - dict2[key]

# padding 矩阵
def Padding(matrix: np.ndarray, wrap:int) -> np.ndarray:
    shape = matrix.shape
    padding = np.zeros((shape[0]+shape[1]-2, shape[0] + shape[1]-2))
    padding[:(shape[0]-1),:(shape[1]-1)] = matrix[:-1,:-1]
    padding[(shape[0]-1):,:(shape[1]-1)] = wrap
    padding[:(shape[0]-1),(shape[1]-1):] = wrap
    for i in range(shape[0]-1):
        padding[i,shape[1]-1+i] = matrix[i,-1]
    for i in range(shape[1]-1):
        padding[shape[0]-1+i,i] = matrix[-1,i]
    return padding

# unpadding 矩阵
def Unpadding(matrix: np.ndarray, dots0: int, dots1: int) -> np.ndarray:
    unpadding = np.zeros((dots0+1, dots1+1))
    unpadding[:-1,:-1] = matrix[:dots0,:dots1]
    for i in range(dots1):
        unpadding[-1,i] = matrix[dots0+i,i]
    for i in range(dots0):
        unpadding[i,-1] = matrix[i,dots1+i]
    return unpadding


# 按统一格式解析命令行参数
def ArgParse():
    global dot_sub, dot_ins, dot_del, edge_sub, edge_ins, edge_del, root_path, void, inf

    parser = argparse.ArgumentParser()
    parser.add_argument("dot_sub", help="Substitution cost of a dot", type=int)
    parser.add_argument("dot_del_ins", help="Deletion/Insertion cost of a dot", type=int)
    parser.add_argument("edge_sub", help="Substitution cost of an edge", type=int)
    parser.add_argument("edge_del_ins", help="Deletion/Insertion cost of an edge", type=int)
    parser.add_argument("graph1", help="Use a graph in Alkane as g_1")
    parser.add_argument("graph2", help="Use a graph in Alkane as g_2")
    parser.add_argument("-m", "--muta", action="store_true",
                        help="Specify whether use Alkane or MUTA. Alkane by default.")
    args = parser.parse_args()

    dot_sub = args.dot_sub
    dot_ins = args.dot_del_ins
    dot_del = args.dot_del_ins
    edge_sub = min(args.edge_sub, 2 * args.edge_del_ins)  # 边替换尝试用边删除+边插入替代
    edge_ins = args.edge_del_ins
    edge_del = args.edge_del_ins

    root_path = root_path + ('MUTA' if args.muta else 'alkane')

    return dot_sub, dot_ins, dot_del, edge_sub, edge_ins, edge_del, root_path, void, inf, \
           args.graph1, args.graph2


class Graph:
    def __init__(self, path: str):
        with open(path, 'r') as fp:
            line = fp.readlines()
            for i in range(len(line)):
                line[i] = line[i].split()
            self.dots, self.edges = tuple(map(int, line[0]))
            self.dtag = [0 for _ in range(self.dots)]
            self.etag = {}
            self.adja = [[] for _ in range(self.dots)]
            self.etagset = set()

            for dot in range(self.dots):
                self.dtag[int(line[dot + 1][0]) - 1] = line[dot + 1][1]
            for i in range(1 + self.dots, 1 + self.dots + self.edges):
                self.adja[int(line[i][0]) - 1].append((int(line[i][1]) - 1, line[i][2]))
                self.adja[int(line[i][1]) - 1].append((int(line[i][0]) - 1, line[i][2]))
                self.etag[(int(line[i][0]) - 1, int(line[i][1]) - 1)] = line[i][2]
                self.etag[(int(line[i][1]) - 1, int(line[i][0]) - 1)] = line[i][2]
            self.etagset = set(self.etag.values())


# 给定点的对应关系, 计算相应的GED
def ComputeGED(graph1: Graph, graph2: Graph, map: tuple) -> int:
    global dot_sub, dot_ins, dot_del, edge_sub, edge_ins, edge_del

    result = 0

    # 点编辑代价
    for dot1 in range(graph1.dots):
        if map[dot1] == void:
            result += dot_del
        else:
            result += (dot_sub if graph1.dtag[dot1] != graph2.dtag[map[dot1]] else 0)
    image = set(map)
    result += dot_ins * len(set(range(graph2.dots)) - image)

    # 边编辑代价
    for dot1 in range(graph1.dots):
        for dot2, tag in graph1.adja[dot1]:
            if dot2 < dot1:  # 每条边只访问一次
                continue
            dot1_map = map[dot1]
            dot2_map = map[dot2]
            if (dot2_map, dot1_map) in graph2.etag:
                # 肯定做边替换，而非边删除+边插入
                result += (edge_sub if graph2.etag[(dot2_map, dot1_map)] != tag else 0)
            else:
                result += edge_del
    for dot1 in range(graph2.dots):
        for dot2, tag in graph2.adja[dot1]:
            if dot2 < dot1:  # 每条边只访问一次
                continue
            if dot1 not in image or dot2 not in image:
                result += edge_ins
            elif (map.index(dot1), map.index(dot2)) not in graph1.etag:
                result += edge_ins
    return result

# 使用 funct 作为点编辑导致的边编辑的代价估计函数
# 产生的是 (n+1)*(m+1) 阶代价矩阵
def ReducedCostMatrix(graph1: Graph, graph2: Graph, funct: 'function') -> np.ndarray:
    global dot_sub, dot_ins, dot_del

    cost_mat = np.zeros((graph1.dots + 1, graph2.dots + 1))
    # 构造左上角和右上角
    for i in range(graph1.dots):
        for j in range(graph2.dots):
            cost_mat[i, j] = (0 if graph1.dtag[i] == graph2.dtag[j] else dot_sub) + funct(graph1, graph2, i, j)
    # 构造右上角和左下角
    cost_mat[-1,:-1] = dot_ins
    cost_mat[:-1, -1] = dot_del

    return cost_mat

# 使用 funct 作为点编辑导致的边编辑的代价估计函数
# 产生的是 (n+m)*(n+m) 阶代价矩阵
def CostMatrix(graph1: Graph, graph2: Graph, funct: 'function') -> np.ndarray:
    global inf
    return Padding(ReducedCostMatrix(graph1, graph2, funct), inf)


# 统一格式打印结果
def PrintGED(graph1: Graph, graph2: Graph, mapping: tuple):
    print("GED =", ComputeGED(graph1, graph2, tuple(mapping)))
    for i in range(graph1.dots):
        print(i + 1, "->", (mapping[i] + 1) if isinstance(mapping[i], int) else None)


# 考虑 1-star 局部结构时额外增加的代价
def Add1Star(graph1: Graph, graph2: Graph, dot1: int, dot2: int) -> int:
    global edge_sub, edge_ins, edge_del

    if dot1 == void and dot2 == void:  # void -> void
        return 0
    elif dot1 == void:  # void -> u
        return edge_ins * len(graph2.adja[dot2])
    elif dot2 == void:  # v -> void
        return edge_del * len(graph1.adja[dot1])
    else:  # v -> u
        edgeset1 = {}
        edgeset2 = {}
        for _, tag in graph1.adja[dot1]:
            edgeset1[tag] = edgeset1.get(tag, 0) + 1
        for _, tag in graph2.adja[dot2]:
            edgeset2[tag] = edgeset2.get(tag, 0) + 1
        InplaceDifferDict(edgeset1, edgeset2)
        # 计算差异集大小
        edge1 = 0
        edge2 = 0
        for val in edgeset1.values():
            if val > 0:
                edge1 += val
            else:
                edge2 -= val
        # 边替换总不劣于边删除+边插入
        if edge1 > edge2:
            return edge2 * edge_sub + (edge1 - edge2) * edge_del
        else:
            return edge1 * edge_sub + (edge2 - edge1) * edge_ins


# 解线性指派问题
def SolveLSAP(cost_matrix: np.ndarray):
    row, col = linear_sum_assignment(cost_matrix)
    # 构造置换矩阵
    permute_mat = np.zeros_like(cost_matrix)
    permute_mat[row, col] = 1
    return row, col, permute_mat


# 构造线性变换的和二次型的矩阵 cost, delta
def ReducedDeltaMatrix(graph1: Graph, graph2: Graph):
    # 压缩的线性变换的矩阵
    cost = ReducedCostMatrix(graph1, graph2, lambda x, y, z, w: 0)
    # 二次型的矩阵
    tot = (graph1.dots + 1) * (graph2.dots + 1)
    delta = np.zeros((tot, tot))
    for i in range(tot):
        for j in range(tot): # dot1 -> dot1_map; dot2 -> dot2_map
            dot1, dot1_map = i // (graph2.dots + 1), i % (graph2.dots + 1)
            dot2, dot2_map = j // (graph2.dots + 1), j % (graph2.dots + 1)
            # 根据对应边的存在性决定边编辑的代价
            if (dot1, dot2) in graph1.etag:
                if (dot1_map, dot2_map) in graph2.etag:
                    delta[i,j] = (edge_sub if graph1.etag[(dot1, dot2)] != graph2.etag[(dot1_map, dot2_map)] else 0)
                else:
                    delta[i,j] = edge_del
            else:
                delta[i,j] = (edge_ins if (dot1_map, dot2_map) in graph2.etag else 0)
    return cost, delta

if __name__ == '__main__':
    dot_ins = 3
    dot_del = 5
    edge_ins = 7
    edge_del = 9
    g1 = Graph(r'C:\Users\lenovo\PycharmProjects\untitled1\GED_data\preprocessed_C\alkane\molecule002.txt')
    g2 = Graph(r'C:\Users\lenovo\PycharmProjects\untitled1\GED_data\preprocessed_C\alkane\molecule003.txt')
    C = ReducedCostMatrix(g1, g2, lambda x,y,z,w: 0)
    print(C)
    print(Padding(C, -1))
    print(Unpadding(Padding(C, -1), g1.dots, g2.dots))
