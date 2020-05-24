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
本次项目中我们进行了诸多连贯的尝试，因此出现了大量的代码复用现象。  
为了使代码的逻辑层次更清楚，以免自己陷入无意义的重复之中，我们对文件进行模块化处理。  
现在，特殊重要的代码都集中在 common 模块中，它们将在若干文件中被反复调用。  
这其中包括表示化学分子的类和 LSAP Solver 等。使用 common.\_\_doc\_\_ 属性查看详细说明。  

## 交互
所有程序均设计为命令行交互，并且由 common 模块的 ArgParse 函数提供了统一的命令行格式。具体格式参见 [ICPR 2016](https://gdc2016.greyc.fr/)。  
为了方便选择测试 Alkane 还是 MUTA，额外提供了可选参数 -m (--muta)。默认将在 Alkane 上测试，指定 -m (--muta) 表明在 MUTA 上测试。  
使用 -h 查看详细说明。  
  
示例：\>\> python BGM1.py -m 2 4 1 1 molecule_100.txt molecule_1000.txt  
&emsp;&emsp;&emsp;\>\> python IFIP.py 2 4 1 2 molecule100.txt molecule010.txt

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
