# -*- coding:utf-8 -*-

"""
BGM1.py: Bipartite Graph Matching with 1-star local structure

Written by Ding Rui
Latest Version: 2020/5/24
"""

from common import ArgParse, Graph, CostMatrix, PrintGED, Add1Star, SolveLSAP

dot_sub, dot_ins, dot_del, edge_sub, edge_ins, edge_del, root_path, void, inf,\
    g1, g2 = ArgParse()


g1 = Graph(root_path + '/'+ g1)
g2 = Graph(root_path + '/' + g2)

_, col, _ = SolveLSAP(CostMatrix(g1, g2, Add1Star))
answer = [(int(col[i]) if col[i] < g2.dots else void) for i in range(g1.dots)]
PrintGED(g1, g2, tuple(answer))
