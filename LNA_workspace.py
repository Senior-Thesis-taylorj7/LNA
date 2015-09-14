import math
import numpy as np
P = [[-2,1],[0,2],[2,1],[0,0],[-1,-2]]

def normDiff(p1,p2):
    return pow(pow(p1[0]-p2[0],2)+pow(p1[1]-p2[1],2),.5)

def gaussian(i,j,p,sigma):
    if abs(i-j) > 1:
        return round(pow(math.e,-(normDiff(p[i-1],p[j-1])*normDiff(p[i-1],p[j-1]))/(sigma*sigma)),2)
    else:
        return 0

def Diagonal(p,sigma,i):
    s = 0
    for j in range(1,len(P)+1):
        s = s + gaussian(i,j,p,sigma)
    return s

def L(p,sigma):
    size = len(p)
    I = np.identity(size)
    Omega = np.zeros([size,size])
    for i in range(1,size+1):
        for j in range(1,size+1):
            Omega[i-1][j-1] = gaussian(i,j,p,sigma)
    diag = np.zeros([size,size])
    for i in range(1,size+1):
        diag[i-1][i-1] = Diagonal(p,sigma,i)
    L = I - np.linalg.inv(diag)*Omega
    return L

def Norm_Sequence(p,sigma):
    Laplacian = L(p,sigma)
    coords = Laplacian * p
    Psigma = []
    for i in range(0,len(coords)):
        Psigma = Psigma+[normDiff(coords[i],[0,0])]
    return Psigma

def P_Sequence(p,sigma_list):
    P_seq = []
    for i in sigma_list:
        P_seq = P_seq + [Norm_Sequence(p,i)]
    return P_seq
