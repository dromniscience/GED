# -*- coding:utf-8 -*-

from scipy.optimize import linear_sum_assignment
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("dot_sub",help="Substitution cost of a dot",type=int)
parser.add_argument("dot_del_ins",help="Deletion/Insertion cost of a dot",type=int)
parser.add_argument("edge_sub",help="Substitution cost of an edge",type=int)
parser.add_argument("edge_del_ins",help="Deletion/Insertion cost of an edge",type=int)
parser.add_argument("graph1",help="Use a graph in Alkane as g_1")
parser.add_argument("graph2",help="Use a graph in Alkane as g_2")
parser.add_argument("-m","--muta", action="store_true",help="Specify whether use Alkane or MUTA. Alkane by default.")
args = parser.parse_args()

dot_sub = args.dot_sub
dot_ins = args.dot_del_ins
dot_del = args.dot_del_ins
edge_sub = args.edge_sub
edge_ins = args.edge_del_ins
edge_del = args.edge_del_ins

root_path = r'./GED_data/preprocessed_C/' + ('MUTA' if args.muta else 'alkane')
void = None
inf = 0x7FFFFFFF # 当然更大也行。这是一个潜在的风险，更好的做法是由输入推断设置上界。

class Graph:
    def __init__(self, path: str):
        with open(root_path + '/' + path, 'r') as fp:
            line = fp.readlines()
            for i in range(len(line)):
                line[i] = line[i].split()
            self.dots, self.edges = tuple(map(int, line[0]))
            self.dtag = [0 for _ in range(self.dots)]
            self.etag = {}
            self.adja = [[] for _ in range(self.dots)]

            for dot in range(self.dots):
                self.dtag[int(line[dot + 1][0]) - 1] = line[dot + 1][1]
            for i in range(1 + self.dots, 1 + self.dots + self.edges):
                self.adja[int(line[i][0]) - 1].append((int(line[i][1]) - 1, line[i][2]))
                self.adja[int(line[i][1]) - 1].append((int(line[i][0]) - 1, line[i][2]))
                self.etag[(int(line[i][0]) - 1, int(line[i][1]) - 1)] = line[i][2]
                self.etag[(int(line[i][1]) - 1, int(line[i][0]) - 1)] = line[i][2]


# 差异集计算
def inplace_differ_dict(dict1 : dict, dict2: dict):
    for key in dict2:
        if key in dict1:
            dict1[key] -= dict2[key]
        else:
            dict1[key] = -dict2[key]

def additional_edge_cost(graph1: Graph, graph2: Graph, dot1: int, dot2: int) -> int:
    if dot1 == void and dot2 == void: # void -> void
        return 0
    elif dot1 == void: # void -> u
        return edge_ins * len(graph2.adja[dot2])
    elif dot2 == void: # v -> void
        return edge_del * len(graph1.adja[dot1])
    else: # v -> u
        edgeset1 = {}
        edgeset2 = {}
        for _,tag in graph1.adja[dot1]:
            edgeset1[tag] = edgeset1.get(tag, 0) + 1
        for _,tag in graph2.adja[dot2]:
            edgeset2[tag] = edgeset2.get(tag, 0) + 1
        inplace_differ_dict(edgeset1, edgeset2)
        # 计算差异集大小
        edge1 = 0
        edge2 = 0
        for val in edgeset1.values():
            if val > 0:
                edge1 += val
            else:
                edge2 -= val
        if edge_sub >= edge_del + edge_ins:
            return edge1 * edge_del + edge2 * edge_ins
        elif edge1 > edge2:
            return edge2 * edge_sub + (edge1 - edge2) * edge_del
        else:
            return edge1 * edge_sub + (edge2 - edge1) * edge_ins

def compute_GED(graph1: Graph, graph2: Graph, map: tuple) -> int:
    result = 0

    # 点编辑代价
    for dot1 in range(g1.dots):
        if map[dot1] == void:
            result += dot_del
        else:
            result += (dot_sub if graph1.dtag[dot1] != graph2.dtag[map[dot1]] else 0)
    image = set(map)
    result += dot_ins * len(set(range(g2.dots)) - image)

    #边编辑代价
    for dot1 in range(g1.dots):
        for dot2,tag in g1.adja[dot1]:
            if dot2 < dot1: # 每条边只访问一次
                continue
            dot1_map = map[dot1]
            dot2_map = map[dot2]
            if (dot2_map, dot1_map) in graph2.etag:
                result += (edge_sub if graph2.etag[(dot2_map, dot1_map)] != tag else 0)
            else:
                result += edge_del
    for dot1 in range(g2.dots):
        for dot2, tag in g2.adja[dot1]:
            if dot2 < dot1: # 每条边只访问一次
                continue
            if dot1 not in image or dot2 not in image:
                result += edge_ins
            elif (map.index(dot1), map.index(dot2)) not in graph1.etag:
                result += edge_ins
    return result


g1 = Graph(args.graph1)
g2 = Graph(args.graph2)


cost_mat = np.zeros((g1.dots + g2.dots, g1.dots + g2.dots))
for i in range(g1.dots):
    for j in range(g2.dots):
        cost_mat[i, j] = (0 if g1.dtag[i] == g2.dtag[j] else dot_sub) + additional_edge_cost(g1, g2, i, j)
    for j in range(g1.dots):
        cost_mat[i, g2.dots + j] = ((dot_del + additional_edge_cost(g1, g2, i, void)) if i == j else inf)
for i in range(g2.dots):
    for j in range(g2.dots):
        cost_mat[i + g1.dots, j] = ((dot_ins + additional_edge_cost(g1, g2, void, j)) if i == j else inf)


row, col = linear_sum_assignment(cost_mat)
answer = [(int(col[i]) if col[i] < g2.dots else void) for i in range(g1.dots)]

print("GED =", compute_GED(g1, g2, tuple(answer)))
for i in range(g1.dots):
    print(i + 1, "->", (answer[i] + 1) if isinstance(answer[i], int) else None)
