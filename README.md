# GED
GED Graphic Edit Distance - 2020算分课程团队项目

## 数据
来源于比赛官网 [ICPR 2016](https://gdc2016.greyc.fr/)。  
已将 .gxl 格式文件转换为 .json 和 .txt 类型。转换后数据存放在文件夹 data 中。

## 损失函数
我们考虑的损失函数的范围与 [ICPR 2016](https://gdc2016.greyc.fr/) 的 Alkane 和 MUTA 数据集相同。  
我们将针对这种代价不受不同标签相似度影响的损失函数在程序中使用特定的优化。

## A*-algorithm
    基于 A* 算法的精确 GED 程序，实现在 A_Star_Exact.py 中。  
交互： 命令行输入, 具体格式见 [ICPR 2016](https://gdc2016.greyc.fr/)  
示例： \>\> python A_Star_Exact.py 2 4 1 1 molecule100.txt molecule50.txt  

## Bipartite Graph Matching
    基于论文 Approximate graph edit distance computation by means of bipartite graph matching(2008) Kaspar Riesen et al.
