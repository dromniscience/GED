# GED
GED Graphic Edit Distance - 2020算分课程团队项目  
项目运行环境：Python 3.x &nbsp; &nbsp; NumPy 最好1.14.1及以上

## 数据
来源于比赛官网 [ICPR 2016](https://gdc2016.greyc.fr/)。  
已将 .gxl 格式文件转换为 .json 和 .txt 类型。转换后数据存放在文件夹 data 中。

## 损失函数
我们考虑的损失函数的范围与 [ICPR 2016](https://gdc2016.greyc.fr/) 的 Alkane 和 MUTA 数据集相同。  
我们将针对这种代价不受不同标签相似度影响的损失函数在程序中使用特定的优化。

## A*-algorithm
    基于 A* 算法的精确 GED 程序，实现在 A_Star.py 中。  
交互： 命令行输入, 具体格式见 [ICPR 2016](https://gdc2016.greyc.fr/)  
示例： \>\> python A_Star.py 2 4 1 1 molecule100.txt molecule050.txt  

## Bipartite Graph Matching
### LSAP + 1-star local structure
    基于论文 Approximate graph edit distance computation by means of bipartite graph matching(2008) Kaspar Riesen et al.
    实现在 BP_Match.py 中。
交互：  命令行输入， \>\> python BP_Match.py -h 查看交互说明  
示例： \>\> python BP_Match.py -m 2 4 1 1 molecule_100.txt molecule_1000.txt

### QAP
&emsp;使用BP_Match(cost, g<sub>1</sub>, g<sub>2</sub>)初始化的 IPFP 方法。 
