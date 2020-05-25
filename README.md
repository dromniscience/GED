# GED
GED Graphic Edit Distance - 2020算分课程团队项目  
项目运行环境：Python 3.x &nbsp; &nbsp; NumPy 最好1.14.1及以上

## 数据
来源于比赛官网 [ICPR 2016](https://gdc2016.greyc.fr/)。  
已将 .gxl 格式文件转换为 .json 和 .txt 类型。转换后数据存放在文件夹 data 中。  
详细的数据说明请见该目录下 README.md 文件。

## 损失函数
我们考虑的损失函数的范围与 [ICPR 2016](https://gdc2016.greyc.fr/) 的 Alkane 和 MUTA 数据集相同。  
我们将针对这种**代价不受不同标签相似度影响**的损失函数在程序中使用特定的优化。

## 代码组织逻辑
&emsp;&emsp;本次项目中我们进行了诸多连贯的尝试，因此出现了大量的代码复用现象。为了使代码的逻辑层次更清楚，以免自己陷入无意义的重复之中，我们对文件进行模块化处理。  
&emsp;&emsp;现在，特殊重要的代码都集中在 common 模块中，它们将在若干文件中被反复调用。这其中包括表示化学分子的 Graph 类等。使用 common.\_\_doc\_\_ 属性查看详细说明。  

## 交互
&emsp;&emsp;所有程序均设计为命令行交互，并且由 common 模块的 ArgParse 函数提供了统一的命令行格式。具体格式参见 [ICPR 2016](https://gdc2016.greyc.fr/)。  
&emsp;&emsp;为了方便选择测试 Alkane 还是 MUTA，额外提供了可选参数 -m (--muta)。默认将在 Alkane 上测试，指定 -m (--muta) 表明在 MUTA 上测试。使用 -h 查看详细说明。  
  
示例：\>\> python BGM1.py -m 2 4 1 1 molecule_100.txt molecule_1000.txt  
&emsp;&emsp;&emsp;\>\> python IPFP1.py 2 4 1 2 molecule100.txt molecule010.txt

## A*-algorithm
    基于 A* 算法的精确 GED 程序，实现在 AStar.py 中. 

## Bipartite Graph Matching
### LSAP + 1-star local structure
    基于论文 Approximate graph edit distance computation by means of bipartite graph matching(2008) Kaspar Riesen et al.
    实现在 BGM1.py 中。
    
### QAP
    IPFP 算法部分基于论文 A Quadratic Assignment Formulation of the Graph Edit Distance(2015) Sebastien Bougleux et al.
    实现在 IPFP1.py 和 IPFP2.py 中。两个文件的区别仅在 IPFP 算法初始化方法不同。
    此处修正了一个原文伪码的疑似错误: 最后修正 IPFP 算法的答案为 binary bi-stochastic matrix 时，最相近的矩阵应该 argmax 目标值而非 argmin。 

## Addendum
&emsp;&emsp;概括地说，常规的 LSAP 问题的有效解法包括以下三种:
* 最大流
* 最短增广路径
* 线性规划

&emsp;&emsp;尽管原始-对偶规划 (primal-dual algorithm) 常常在实践中表现出最明显的优势，但我们仍然很时有兴趣回到这个问题最初被有效解决的情况。那是得益于 Konig 提供的有趣命题，即“有限简单二部图的最大匹配数恰好等于最小点覆盖数”，所提出的解法。尽管它有着 O(n<sup>4</sup>) 的复杂度，但本身这个算法的设计却非常有趣，还很有娱乐性。它是被 J. Munkres 严谨地描述并证明的 _(Algorithms for the assignment and transportation problems 1957)_ 。那时计算机并不普及，算法的实现常常是手动模拟，因此生动的描述是它能够广泛传播的重要优势。例如原来的算法就是用对零标星 (starred zero) 或者加撇 (primed zero) 讲述的。我们在这里完整地复现了这个有趣的设计。当然，写出来的代码也很美观。  
&emsp;&emsp;代码在 lsap 模块中。当图的顶点不多时，我们会调用自己的 LSAP Solver；当然，它在性能上和 SciPy 的 \_lsap\_module 相比还有明显的差距。
