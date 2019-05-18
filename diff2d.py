import numpy as numpy

class Diffusion_2d(object):

    def __init__(self, t, D=1, Nt=100, Nx=100, Ny=100, x_bound=(0, 1), y_bound=(0, 1)):
        """
        初始化的时候已经确定了网格和模拟的时间.
        """
        # xy的边界. 默认2比1大.
        self.x1, self.x2 = x_bound
        self.y1, self.y2 = y_bound

        # 时间长度
        self.t = t
        # t网格化的点数, 即迭代次数
        self.Nt = Nt

        # xy网格化的点数
        self.Nx = Nx
        self.Ny = Ny

        # 所有的步长
        self.hx = (self.x2 - self.x1)/float(self.Nx - 1)
        self.hy = (self.y2 - self.y1)/float(self.Ny - 1)
        self.ht = t/float(self.Nt-1)

        # 扩散系数
        self.D = D


    def explicit(self, u=np.zeros([Nx, Ny])):
        unew = np.zeros([Nx, Ny])



        