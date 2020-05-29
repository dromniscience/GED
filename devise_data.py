"""
devise_data.py: make devise data for testing

using class Graph in common.py

Written by Jiang Heng
Last Version: 2020/5/29

"""


import numpy as np
import random

# choose cost, finally 2,4,1,2
cost_list = [[2, 4, 1, 1], [2, 4, 1, 2], [6, 2, 2, 1]]
tmp_list = cost_list[1]
dot_sub = tmp_list[0]
dot_ins = tmp_list[1]
dot_del = tmp_list[1]
edge_sub = tmp_list[2]
edge_ins = tmp_list[3]
edge_del = tmp_list[3]

root_path = r'./GED_data/preprocessed_C/MUTA'

out_ans = './dev_data2412/ans.txt'
ans_f = open(out_ans, 'w')


# naive loop for all MUTA documents
for names in range(1, 4338):
    path = 'molecule_' + str(names) + '.txt'
    origin = 'molecule_0.txt'

    g1 = Graph(path)
    g2 = Graph(origin)
    cost = 0

    '''random dots deletion'''

    to_del_dots = random.randint(0, int(g1.dots / 2))
    cost = 0
    extra_edge_del = 0

    # auxiliary list for the dots mapping from g1 to g2
    aux_list = [i for i in range(g1.dots)]
    del_list = random.sample(aux_list, to_del_dots)
    after_list = list(set(aux_list).difference(del_list))

    # relative edge deletion
    for edges in g1.etag:
        if (edges[0] in del_list) or (edges[1] in del_list):
            g1.etag[edges] = "None"
            extra_edge_del += 1

    # make g2.dtag
    for i in after_list:
        g2.dtag.append(g1.dtag[i])

    # make g2.etag(edges not deleted)
    for i in g1.etag.keys():
        if g1.etag[i] != "None":
            g2.etag[(after_list.index(i[0]), after_list.index(i[1]))] = g1.etag[i]

    # make g2 fully constructed(after node deletion)
    g2.dots = g1.dots - to_del_dots
    g2.edges = g2.etag.__len__()

    # calculate cost
    cost += to_del_dots * dot_del
    cost += int(extra_edge_del/2) * edge_del

    '''random edge deletion'''

    to_del_edges = 0
    for i in range(int(g2.edges/2)):
        s = random.randint(0, g2.dots-2)
        t = random.randint(s+1, g2.dots-1)
        if (s, t) in g2.etag:
            del g2.etag[(s, t)]
            del g2.etag[(t, s)]
            to_del_edges += 1

    cost += to_del_edges * edge_del

    '''random substitution'''

    to_sub_dots = random.randint(0, g2.dots)
    to_sub_edges = random.randint(0, int(g2.edges/2))

    # in order to deal with the repeated sub easily.jpg
    real_sub_dots = to_sub_dots
    real_sub_edges = to_sub_edges

    # naive random dots substitution
    for i in range(to_sub_dots):
        pos = random.randint(0, g2.dots - 1)
        if g2.dtag[pos] == "none":
            real_sub_dots -= 1
        else:
            g2.dtag[pos] = "none"

    # naive random edge substitution
    for i in range(to_sub_edges):
        pos = random.choice(list(g2.etag))
        if g2.etag[pos] == "none":
            real_sub_edges -= 1
        else:
            g2.etag[pos] = "none"
            tmp = (pos[1], pos[0])
            g2.etag[tmp] = "none"

    cost += real_sub_dots * dot_sub
    cost += real_sub_edges * edge_sub

    # output
    outpath = './dev_data2412/dev_' + str(names) + '.txt'

    f = open(outpath, "w")

    print("%d %d" % (g2.dots, int(g2.edges)/2), file=f)
    for i in range(g2.dots):
        print("%d %s" % (i+1, g2.dtag[i]), file=f)
    for pair in g2.etag.keys():
        if pair[0] < pair[1]:
            print("%d %d %s" % (pair[0]+1, pair[1]+1, g2.etag[pair]), file=f)
    print(names, ":", cost, file=ans_f)
    f.close()
