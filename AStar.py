# -*- coding:utf-8 -*- 

"""
AStar.py: A*-algorithm for GED

Written by Ding Rui
Last Version: 2020/5/24
"""

from queue import PriorityQueue
from common import ArgParse, Graph, PrintGED, InplaceDifferDict

dot_sub, dot_ins, dot_del, edge_sub, edge_ins, edge_del, root_path, void, inf,\
    graph1, graph2 = ArgParse()

# 局部路径
class Partial:
    def __init__(self, graph1: Graph, graph2: Graph, cost: int, old_part_map=(), expand=void):
        self.graph1, self.graph2 = graph1, graph2
        self.part_map = old_part_map + (expand,)
        self.cost = cost + self.Add_Cost()
        self.tot_cost = self.cost + self.Heuristic()

    def __lt__(self, other):
        return self.tot_cost < other.tot_cost

    def Add_Cost(self) -> int:
        dot1 = len(self.part_map) - 1
        dot2 = self.part_map[-1]
        result = 0

        # 最后加入的是替换还是删除
        if dot2 == void:
            # 点删除
            result += dot_del
            # 边删除(但是只考虑序号更小的边)
            for dot in range(dot1):
                if (dot1, dot) in self.graph1.etag:
                    result += edge_del
        else:
            # 点替换
            result += (dot_sub if self.graph1.dtag[dot1] != self.graph2.dtag[dot2] else 0)
            # 边操作(但是只考虑序号更小的边)
            for dot in range(dot1):
                dot_map = self.part_map[dot]  # dot1 -> dot2; dot -> dot_map;
                if (dot1, dot) in self.graph1.etag:
                    if (dot2, dot_map) in self.graph2.etag:
                        # 肯定做的是边替换(而不是边删除+边插入)
                        result += (edge_sub if self.graph2.etag[(dot2, dot_map)] != self.graph1.etag[(dot1, dot)] else 0)
                    else:
                        result += edge_del
                elif (dot_map, dot2) in self.graph2.etag:
                    result += edge_ins

        # 当前是否能确定graph2剩余点的映射
        if dot1 == self.graph1.dots - 1:
            # 点代价
            already = set(self.part_map) - {void}
            left = set(range(self.graph2.dots)) - already
            result += len(left) * dot_ins
            # 边代价
            for i in left:
                for j in already:
                    if (i,j) in self.graph2.etag:
                        result += edge_ins
                already.add(i)

        return result

    def Heuristic(self) -> int:
        result = 0
        # 余下的点编辑代价的下界
        dot = len(self.part_map)
        # 是否需要启发
        if dot == self.graph1.dots:
            return 0
        dict1 = {tag:0 for tag in self.graph1.dtag[dot:]}
        for tag in self.graph1.dtag[dot:]:
            dict1[tag] += 1

        dotset = set(range(self.graph2.dots)) - set(self.part_map)
        dict2 = {self.graph2.dtag[num]:0 for num in dotset}
        for num in dotset:
            dict2[self.graph2.dtag[num]] += 1
        InplaceDifferDict(dict1, dict2)
        # 获取标签分布的差异集大小
        dots1 = 0
        dots2 = 0
        for val in dict1.values():
            if val > 0:
                dots1 += val
            else:
                dots2 -= val
        # 删除/插入 or 替换
        if dot_sub >= dot_del + dot_ins:
            result += (dots1 * dot_del + dots2 * dot_ins)
        elif dots1 > dots2:
            result += (dots2 * dot_sub + (dots1 - dots2) * dot_del)
        else:
            result += (dots1 * dot_sub + (dots2 - dots1) * dot_ins)

        # 余下的边替换编辑代价的下界
        dict1 = {tag: 0 for tag in self.graph1.etagset}
        dict2 = {tag: 0 for tag in self.graph2.etagset}
        for (st, ed), tag in self.graph1.etag.items():
            if st > ed: # 去掉平行边
                continue
            if st >= dot or ed >= dot:
                dict1[tag] += 1
        for (st, ed), tag in self.graph2.etag.items():
            if st > ed: # 去掉平行边
                continue
            if st in dotset or ed in dotset:
                dict2[tag] += 1
        InplaceDifferDict(dict1, dict2)
        # 获取标签分布的差异集大小
        edges1 = 0
        edges2 = 0
        for val in dict1.values():
            if val > 0:
                edges1 += val
            else:
                edges2 -= val
        # 边替换总不劣于边删除+边插入
        if edges1 > edges2:
            result += (edges2 * edge_sub + (edges1 - edges2) * edge_del)
        else:
            result += (edges1 * edge_sub + (edges2 - edges1) * edge_ins)
        return result


graph1 = Graph(root_path + '/'+ graph1)
graph2 = Graph(root_path + '/'+ graph2)

queue = PriorityQueue()
for i in range(graph2.dots):
    queue.put(Partial(graph1, graph2, 0, tuple(), i))
queue.put(Partial(graph1, graph2, 0, tuple(), void))

while True:
    partial = queue.get()
    if len(partial.part_map) == graph1.dots:
        PrintGED(graph1, graph2, partial.part_map)
        break
    else:
        left = set(range(graph2.dots)) - set(partial.part_map)
        for i in left:
            queue.put(Partial(graph1, graph2, partial.cost, partial.part_map, i))
        queue.put(Partial(graph1, graph2, partial.cost, partial.part_map, void))
