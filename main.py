import numpy as np
import time
from kelvin_helmholtz.tools import Kelvin_Helmholtz
import matplotlib.pyplot as plt


def gridmobile(xspeedmatrix, yspeedmatrix, length, parametermatrix, deltat):
    (a, b) = np.shape(xspeedmatrix)
    for i in range(a):
        for j in range(b):
            parameterbefore = parametermatrix[i, j]
            xspeed = xspeedmatrix[i, j]
            yspeed = yspeedmatrix[i, j]
            deltax = xspeed*deltat
            deltay = yspeed*deltat
            x = int((i+deltax//length) % a)
            y = int((j+deltay//length) % b)
            S1 = (length-deltax % length)*(length-deltay % length)
            S2 = deltax % length*(length-deltay % length)
            S3 = (length-deltax % length)*deltay % length
            S4 = deltax % length*deltay % length
            newparametermatrix = np.zeros((a, b))
            newparametermatrix[x, y] += parameterbefore*S1
            newparametermatrix[(x+1) % a, y] += parameterbefore*S2
            newparametermatrix[x, (y+1) % b] += parameterbefore*S3
            newparametermatrix[(x+1) % a, (y+1) % b] += parameterbefore*S4
    return (newparametermatrix)


def avgr(xspeedmatrix, yspeedmatrix, length, rmatrix, Nmatrix, deltat):
    (a, b) = np.shape(xspeedmatrix)
    for i in range(a):
        for j in range(b):
            rbefore = rmatrix[i, j]
            Nbefore = Nmatrix[i, j]
            xspeed = xspeedmatrix[i, j]
            yspeed = yspeedmatrix[i, j]
            deltax = xspeed*deltat
            deltay = yspeed*deltat
            x = int((i+deltax//length) % a)
            y = int((j+deltay//length) % b)
            S1 = (1-deltax % length)*(1-deltay % length)
            S2 = deltax % length*(1-deltay % length)
            S3 = (1-deltax % length)*deltay % length
            S4 = deltax % length*deltay % length
            newrmatrix = np.zeros(x, y)
            newNmatrix = np.zeros(x, y)
            newrmatrix[x, y] += rbefore**3*S1*Nbefore
            newrmatrix[(x+1) % a, j] += rbefore**3*S2*Nbefore
            newrmatrix[i, (j+1) % b] += rbefore**3*S3*Nbefore
            newrmatrix[(i+1) % a, (j+1) % b] += rbefore**3*S4*Nbefore
            newNmatrix[i, j] += Nbefore*S1
            newNmatrix[(i+1) % a, j] += Nbefore*S2
            newNmatrix[i, (j+1) % b] += Nbefore*S3
            newNmatrix[(i+1) % a, (j+1) % b] += Nbefore*S4
    for i in range(x):
        for j in range(y):
            newrmatrix[i][j] = (newrmatrix[i][j]/newNmatrix[i][j])**(1/3)
    return (newrmatrix, newNmatrix)

def main():
    K = Kelvin_Helmholtz()
    K.P = 0*K.X
    t = 0
    for i in range(10):
        print("t={}".format(t))
        t += K.dt
        K.P=gridmobile(K.vx,K.vy,1,K.P,K.dt)
    
        plt.clf()
        plt.imshow(K.P.T)
        plt.clim(0.8, 2.2)
        ax = plt.gca()
        ax.invert_yaxis()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        # plt.draw()
        # time.sleep(0.01)
        plt.pause(0.01)
    plt.show()

print("over")
    # print(K.P)

if __name__ == "__main__":
    main()