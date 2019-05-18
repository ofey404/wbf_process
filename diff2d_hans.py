import numpy as np
import scipy.linalg

def solver_dense(
    I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5, user_action=None):
    
    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)

    dx = x[1] - x[0]
    dy = y[1] - y[0]

    dt = float(dt)
    Nt = int(round(T/dt))

    t = np.linspace(0, Nt*dt, Nt+1)

    Fx = a*dt/dx**2
    Fy = a*dt/dy**2

    u = np.zeros((Nx+1, Ny+1))  # 初始化下一个时间段的u
    u_n = np.zeros((Nx+1, Ny+1))  # 已知的时间段的u.

    Ix = range(0, Nx+1)
    Iy = range(0, Ny+1)
    It = range(0, Nt+1)

    for i in Ix:
        for j in Iy:
            u_n[i, j] = I(x[i], y[j])

    # filling of A
    # A 和时间无关, 是一个常量, 所以只需要初始化一次.

    N = (Nx+1)*(Ny+1)
    A = np.zeros((N, N))
    b = np.zeros(N)

    # 这是什么意思?
    m = lambda i, j: j*(Nx+1) + i

    j = 0
    for i in Ix:
        p = m(i, j)
        A[p, p] = 1

    for j in Iy[1: -1]:
        # 边界
        i = 0
        p = m(i, j)
        A[p, p] = 1
        for i in Ix[1: -1]:
            p = m(i, j)
            A[p, m(i, j-1)] = -theta*Fy
            A[p, m(i-1, j)] = -theta*Fx
            A[p, p] = 1 + 2*theta*(Fx+Fy)
            A[p, m(i+1, j)] = -theta*Fx
            A[p, m(i, j+1)] = -theta*Fy
        i = Nx
        p = m(i, j)
        A[p, p] = 1   # 边界条件

    j = Ny
    for i in Ix:
        p = m(i, j)
        A[p, p] = 1


    ## 利用A迭代.

    for n in It[0:-1]:
        # 计算b
        j = 0
        for i in Ix:
            # 边界条件
            p = m(i, j)
            b[p] = 0
        for j in Iy[1:-1]:
            # 边界条件不变
            i = 0
            p = m(i, j)
            b[p] = 0
            for i in Ix[1:-1]:
                p = m(i, j)
                b[p] = u_n[i, j] + \
                (1-theta) * (
                    Fx*(u_n[i+1, j] - 2*u_n[i, j] + u_n[i-1, j]) +\
                    Fy*(u_n[i, j+1] - 2*u_n[i, j] + u_n[i, j-1])
                ) + \
                theta*dt*f(i*dx, j*dy, (n+1)*dt) + \
                (1 - theta)*dt*f(i*dx, j*dy, n*dt)
            # 边界
            i = Nx
            p = m(i, j)
            b[p] = 0
        j = Ny
        for i in Ix:
            # 边界
            p = m(i, j)
            b[p] = 0

        c = scipy.linalg.solve(A, b)

        # fill u with vector c
        for i in Ix:
            for j in Iy:
                u[i, j] = c[m(i, j)]
                

        u_n, u = u, u_n

    return u_n




