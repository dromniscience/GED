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
Version 1.0: 2020/5/24
Version 1.1: 2020/6/2   Speed up the formulation of the quadratic matrix
'''

# from Common import * 时仅允许导入以下变量
# 但建议不要使用此格式导入本模块!
__all__ = ['ArgParse', 'Graph', 'ComputeGED', 'CostMatrix',
           'PrintGED', 'InplaceDifferDict', 'Add1Star', 'SolveLSAP', 'DeltaMatrix']

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
def CostMatrix(graph1: Graph, graph2: Graph, funct: 'function') -> np.ndarray:
    global dot_sub, dot_ins, dot_del

    cost_mat = np.zeros((graph1.dots + graph2.dots, graph1.dots + graph2.dots))
    # 构造左上角和右上角
    for i in range(graph1.dots):
        for j in range(graph2.dots):
            cost_mat[i, j] = (0 if graph1.dtag[i] == graph2.dtag[j] else dot_sub) + funct(graph1, graph2, i, j)
        for j in range(graph1.dots):
            cost_mat[i, graph2.dots + j] = ((dot_del + funct(graph1, graph2, i, void)) if i == j else inf)
    # 构造左下角
    for i in range(graph2.dots):
        for j in range(graph2.dots):
            cost_mat[i + graph1.dots, j] = ((dot_ins + funct(graph1, graph2, void, j)) if i == j else inf)
    return cost_mat


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
def DeltaMatrix(graph1: Graph, graph2: Graph):
    # 压缩的线性变换的矩阵
    cost = CostMatrix(graph1, graph2, lambda x, y, z, w: 0)
    # 二次型的矩阵
    tot = graph1.dots + graph2.dots
    delta = np.zeros((tot**2, tot**2))
    # 填充
    for dot1, dot2 in graph1.etag:
        delta[dot1*tot:(dot1+1)*tot,dot2*tot:(dot2+1)*tot] = edge_del
    for dot1, dot2 in graph2.etag:
        delta[dot1::tot,dot2::tot] = edge_ins
    for (dot1, dot2), tag1 in graph1.etag.items():
        for (dot1_map, dot2_map), tag2 in graph2.etag.items():
            delta[dot1*tot+dot1_map, dot2*tot+dot2_map] = (edge_sub if tag1 != tag2 else 0)
    for i in range(tot):
        delta[np.diag_indices(tot**2)] = 0
    # 优先级最高, 故最后填充
    for i in range(graph1.dots):
        for j in range(graph1.dots):
            if i == j:
                continue
            delta[i*tot+j+graph2.dots,:] = inf
            delta[:,i*tot+j+graph2.dots] = inf
    for i in range(graph2.dots):
        for j in range(graph2.dots):
            if i == j:
                continue
            delta[(i+graph1.dots)*tot+j,:] = inf
            delta[:,(i+graph1.dots)*tot+j] = inf
    return cost, delta
