import numpy as np
def gridmobile(xspeedmatrix,yspeedmatrix,length,parametermatrix,deltat):
    (x,y)=np.shape(xspeedmatrix)
    for i in range(x):
        for j in range(y):
            parameterbefore=parametermatrix(i,j)
            xspeed=xspeedmatrix[i,j]
            yspeed=yspeedmatrix[i,j]
            deltax=xspeed*deltat
            deltay=yspeed*deltat
            x=i+deltax//length
            y=j+deltay//length
            S1=(1-deltax%length)*(1-deltay%length)
            S2=deltax%length*(1-deltay%length)
            S3=(1-deltax%length)*deltay%length
            S4=deltax%length*deltay%length
            newparametermatrix= np.zeros(x,y)
            newparametermatrix[i,j]+=parameterbefore*S1
            newparametermatrix[i+1,j]+=parameterbefore*S2
            newparametermatrix[i,j+1]+=parameterbefore*S3
            newparametermatrix[i+1,j+1]+=parameterbefore*S4
    return (newparametermatrix)
def avgr(xspeedmatrix,yspeedmatrix,length,rmatrix,Nmatrix,deltat):
    (x,y)=np.shape(xspeedmatrix)
    for i in range(x):
        for j in range(y):
            rbefore=rmatrix[i,j]
            Nbefore=Nmatrix[i,j]
            xspeed=xspeedmatrix[i,j]
            yspeed=yspeedmatrix[i,j]
            deltax=xspeed*deltat
            deltay=yspeed*deltat
            x=i+deltax//length
            y=j+deltay//length
            S1=(1-deltax%length)*(1-deltay%length)
            S2=deltax%length*(1-deltay%length)
            S3=(1-deltax%length)*deltay%length
            S4=deltax%length*deltay%length
            newrmatrix= np.zeros(x,y)
            newNmatrix=np.zeros(x,y)
            newrmatrix[i,j]+=rbefore**3*S1*Nbefore
            newrmatrix[i+1,j]+=rbefore**3*S2*Nbefore
            newrmatrix[i,j+1]+=rbefore**3*S3*Nbefore
            newrmatrix[i+1,j+1]+=rbefore**3*S4*Nbefore
            newNmatrix[i,j]+=Nbefore*S1
            newNmatrix[i+1,j]+=Nbefore*S2
            newNmatrix[i,j+1]+=Nbefore*S3
            newNmatrix[i+1,j+1]+=Nbefore*S4
    for i in range(x):
        for j in range(y):
            newrmatrix[i][j]=(newrmatrix[i][j]/newNmatrix[i][j])**(1/3)
    return (newrmatrix,newNmatrix)

