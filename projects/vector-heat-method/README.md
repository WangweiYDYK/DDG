# Heat Method and Vector Heat Method

## heat method

$\phi(x, y)=\lim _{t \rightarrow 0} \sqrt{-4 t \log k_{t}(x, y)}$

![Alt text](<P8OIJE8_%MQE`((6ALTAAWX.png>)
当热核函数存在误差的时候，直接使用Varadhan公式会有非常显著的误差。

### Alg
在曲面上给定一点x上施加一个热源，并对其进行扩散获得温度场

1.$\frac{d}{dt}u = \Delta u$
取u为顶点坐标的函数，代入热传导方程并离散化后有$\frac{u^{k+1} - u^k}{h} = Lu^k$,隐式迭代下为$(I-hL)u^{k+1} = u^k$

2.$X = -\frac{\nabla u}{|\nabla u|}$
对温度场的梯度进行归一化后得到距离场扩散的向量场X

3.求解$\Delta u = \nabla \cdot X$
对于已知向量场X，如果希望寻找一个势场u使得u可以表示X，则可以构造迪利克雷能量
$E(u)=\int_{M}|\nabla u-X|^{2} d A$
其等价于求解泊松方程$\Delta u = \nabla \cdot X$

### cotLaplacain推导
(https://zhuanlan.zhihu.com/p/372670140)

### 处理点云模型或者多边形模型

### 边界条件处理
1.纽曼边界
2.迪利克雷边界

### 收敛性


### heat method的并行

从热扩散方程中获得一个接近于单位向量场的可积分的梯度场，对该梯度场进行积分获得测地线距离

Gauss-Seidel热扩散和测地线距离积分


### 可优化的部分
时间t的选取
加权距离计算

## Vector Heat Method

### Vector Heat Method的并行分析

### log map

### Karcher Means and Geometric Medians