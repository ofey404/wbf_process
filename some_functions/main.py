import numpy as np
def gridmobile(xspeedmatrix,yspeedmatrix,length,parametermatrix,deltat):
    (a,b)=np.shape(xspeedmatrix)
    for i in range(a):
        for j in range(b):
            parameterbefore=parametermatrix(i,j)
            xspeed=xspeedmatrix[i,j]
            yspeed=yspeedmatrix[i,j]
            deltax=xspeed*deltat
            deltay=yspeed*deltat
            x=(i+deltax//length)%a
            y=(j+deltay//length)%b
            S1=(length-deltax%length)*(length-deltay%length)
            S2=deltax%length*(length-deltay%length)
            S3=(length-deltax%length)*deltay%length
            S4=deltax%length*deltay%length
            newparametermatrix= np.zeros(x,y)
            newparametermatrix[x,y]+=parameterbefore*S1
            newparametermatrix[(x+1)%a,y]+=parameterbefore*S2
            newparametermatrix[x,(y+1)%b]+=parameterbefore*S3
            newparametermatrix[(x+1)%a,(y+1)%b]+=parameterbefore*S4
    return (newparametermatrix)
def avgr(xspeedmatrix,yspeedmatrix,length,rmatrix,Nmatrix,deltat):
    (a,b)=np.shape(xspeedmatrix)
    for i in range(a):
        for j in range(b):
            rbefore=rmatrix[i,j]
            Nbefore=Nmatrix[i,j]
            xspeed=xspeedmatrix[i,j]
            yspeed=yspeedmatrix[i,j]
            deltax=xspeed*deltat
            deltay=yspeed*deltat
            x=(i+deltax//length)%a
            y=(j+deltay//length)%b
            S1=(1-deltax%length)*(1-deltay%length)
            S2=deltax%length*(1-deltay%length)
            S3=(1-deltax%length)*deltay%length
            S4=deltax%length*deltay%length
            newrmatrix= np.zeros(x,y)
            newNmatrix=np.zeros(x,y)
            newrmatrix[x,y]+=rbefore**3*S1*Nbefore
            newrmatrix[(x+1)%a,j]+=rbefore**3*S2*Nbefore
            newrmatrix[i,(j+1)%b]+=rbefore**3*S3*Nbefore
            newrmatrix[(i+1)%a,(j+1)%b]+=rbefore**3*S4*Nbefore
            newNmatrix[i,j]+=Nbefore*S1
            newNmatrix[(i+1)%a,j]+=Nbefore*S2
            newNmatrix[i,(j+1)%b]+=Nbefore*S3
            newNmatrix[(i+1)%a,(j+1)%b]+=Nbefore*S4
    for i in range(x):
        for j in range(y):
            newrmatrix[i][j]=(newrmatrix[i][j]/newNmatrix[i][j])**(1/3)
    return (newrmatrix,newNmatrix)

