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
    L = I - np.dot(np.linalg.inv(diag),Omega)
    return L

def Norm_Sequence(p,sigma):
    Laplacian = L(p,sigma)
    coords = np.dot(Laplacian,p)
    Psigma = []
    for i in range(0,len(coords)):
        Psigma = Psigma+[round(normDiff(coords[i],[0,0]),2)]
    return Psigma

def P_Sequence(p,sigma_list):
    P_seq = []
    for i in sigma_list:
        P_seq = P_seq + [Norm_Sequence(p,i)]
    return P_seq

def S(p,i,q,j,sigma_list,nu):
    if j == (len(q)-1) and i == 0:
        return 0
    if j == 0 and i == (len(p)-1):
        return 0
    if j == 0 and i == 0:
        return 0
    ls = []
    if i > 0:
        ls = ls + [S(p,i-1,q,j,sigma_list)]
    if j > 0:
        ls = ls + [S(p,i,q,j-1,sigma_list)]
    if i > 0 and j > 0:
        ls = ls + [S(p,i-1,q,j-1,sigma_list) + math.pow(math.e,nu*TauGen(p,q,sigma_list)[i][j])]
    return max(ls)

def LNAawk(p,q,sigma_list,nu):
    return S(p,len(p)-1,q,len(q)-1,sigma_list,nu)/math.pow((len(p)-1)*(len(q)-1),.5)

def S_single(p,i,q,j,sigma_list,nu,gap):
    if i == 0 and j == len(q)-1:
        return 0
    if i == 0 and j == 0:
        return 0
    if i == len(p)-1 and j == 0:
        return 0
    ls = [0]
    if i > 0:
        if j > 0:
            ls = ls + [S_single(p,i-1,q,j-1,sigma_list,nu,gap)+1-nu*TauGen(p,q,sigma_list)[i][j]]
        ls = ls + [S_single(p,i-1,q,j,sigma_list,nu,gap)+g]
    if j > 0:
        ls = ls + [S_single(p,i,q,j-1,sigma_list,nu,gap)+g]
    return max(ls)

def LNAswk(p,q,sigma_list,nu,gap):
    return max([S_single(p,i,q,j,sigma_list,nu,gap) for i in range(0,len(p)) for j in range(0,len(q))])