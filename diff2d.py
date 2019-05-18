import numpy as np

class Diffusion_2d(object):

    def __init__(self, Nx=100, Ny=100, x_bound=(0, 1), y_bound=(0, 1)):
        """
        初始化的时候已经确定了
        """
        # xy的边界. 默认2比1大.

        # xy网格化的点数
        self.Nx = Nx
        self.Ny = Ny

        # 所有的步长
        self.hx = (self.x2 - self.x1)/(self.Nx - 1)
        self.hy = (self.y2 - self.y1)/(self.Ny - 1)


    def explicit(self, y_bound=(0, 1), x_bound=(0, 1), u=np.zeros([100, 100]), alpha=1, t=10, Nt=100):
        """Args:
          u: 初始化的密度矩阵. 四边认为是边界条件. 传入的矩阵形状就认为是网格的画法.
          alpha: 扩散方程中的系数
          t: 模拟时长
          Nt: 模拟的轮数
        Returns:
          u_final: t时刻的密度矩阵.
        """
        y1, y2 = y_bound
        x1, x2 = x_bound
        ht = t/(Nt-1)  # t的步长
        Ny, Nx = u.shape

        hx = (x2 - x1)/(Nx - 1)
        hy = (y2 - y1)/(Ny - 1)

        F = 
        

        # 这个条件是hemmar大哥代码里加的, 还没查到文献.
        if ht > 1/2 * (((hx*hy)**2)/(hx**2+hy**2)):
            print("The numerical scheme is not stable for this condition")

        unew = np.zeros([Nx, Ny])

        # 定义一些slice对象, 让代码更可读
        f2tl = slice(None, -2)  # "first to thirdlast" 第一个到第n-2个
        s2sl = slice(1,-1)  # "second to secondlast" 第二个到第n-1个数, 从一开始数.slice
        t2l = slice(2, None)  # "third to last" 第三个到第n个



        for t in range(1, Nt+1):
            unew[s2sl, s2sl] = u[s2sl, s2sl] + 


        
        



        