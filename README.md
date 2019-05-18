

## 一维扩散方程的求解

[The 1D diffusion equation](https://hplgit.github.io/num-methods-for-PDEs/doc/pub/diffu/sphinx/._main_diffu001.html)

时间上需要一个边界条件, 空间上需要两个.

### Forward euler scheme

离散化

把求导换成有限差分. $u_i^n$对应位置为i, 时间分块为n的分量.

$$
u_{i}^{n+1}=u_{i}^{n}+F\left(u_{i+1}^{n}-2 u_{i}^{n}+u_{i-1}^{n}\right)
$$

算法: 知道$u_{any_i}^n$, 计算出$u_i^{n+1}$.

### backward eular scheme

假设只有0123四个点

$$
(1+2 F) u_{1}^{n}-F u_{2}^{n}=u_{1}^{n-1}
$$

$$
-F u_{1}^{n}+(1+2 F) u_{2}^{n}=u_{2}^{n-1}
$$

有N个点就有N-2个方程. 求解联立方程. 隐式欧拉方法更加全面, 它在每个时间点求解耦合方程.

一维联立方程的迭代是不断左乘一个稀疏矩阵.

### Sparse matrix implementation

scipy.sparse包可以解稀疏矩阵的问题.

### Crank-Nicolson scheme

$$
u_{i}^{n+1}-\frac{1}{2} F\left(u_{i-1}^{n+1}-2 u_{i}^{n+1}+u_{i+1}^{n+1}\right)=u_{i}^{n}+\frac{1}{2} F\left(u_{i-1}^{n}-2 u_{i}^{n}+u_{i+1}^{n}\right)
$$

又一种算法.