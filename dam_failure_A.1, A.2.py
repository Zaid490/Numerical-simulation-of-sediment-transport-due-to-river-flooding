# This program solves 1-D Saint-Venant-Exner equations in paper 'Numerical madelling of solid transport caused by an extreme flood: case of the Talal Dam failure

from math import*

from tkinter import *
import copy

#Some Matrix Operations
#matrix multipliction
def matrix_mul2(X,V):#,  mod): #V: vector, X:matrix      ---- not checked!
    L=len(V)
    res=[0]*L
    for i in range(L):
        a=X[i]
        for j in range(L):
            res[i]+=(a[j]*V[j])#%mod
    return res


#inverse of matrix
def eliminate(r1, r2, col, target=0):
    fac = (r2[col]-target) / r1[col]
    for i in range(len(r2)):
        r2[i] -= fac * r1[i]

def gauss(a):
    for i in range(len(a)):
        if a[i][i] == 0:
            for j in range(i+1, len(a)):
                if a[i][j] != 0:
                    a[i], a[j] = a[j], a[i]
                    break
            else:
                raise ValueError("Matrix is not invertible")
        for j in range(i+1, len(a)):
            eliminate(a[i], a[j], i)
    for i in range(len(a)-1, -1, -1):
        for j in range(i-1, -1, -1):
            eliminate(a[i], a[j], i)
    for i in range(len(a)):
        eliminate(a[i], a[i], i, target=1)
    return a

def stabilize(m):
    M=[]
    for j0 in m:
        for j1 in range(len(j0)):
            if j0[j1]<10**-12:
                j0[j1]=0
        M.append(j0)
    return M

def inverse(a):
    tmp = [[] for _ in a]
    for i,row in enumerate(a):
        assert len(row) == len(a)
        tmp[i].extend(row + [0]*i + [1] + [0]*(len(a)-i-1))
    gauss(tmp)
    ret = []
    for i in range(len(tmp)):
        ret.append(tmp[i][len(tmp[i])//2:])
    return stabilize(ret)




#Parameters in formulas A.14, A.15
def B(i, n): #from the regression given in Figure-4
    return 0.0019 * h[i]**3 + 0.0094 * h[i]**2 + 1.7163 * h[i] 

def CG(i, n):
    return teta

def CI(i, n):
    return teta

def CH(i, n):
    return  delx / (4*delt) * ( B(i+1, n) + B(i, n) )

def CJ(i, n):
    return -delx / (4*delt) * ( B(i+1, n) + B(i, n) )

def CK(i, n):
    return -( Q[i+1, n] - Q[i, n] )


def S(i,n): #!
    return w[i]


def Rh(i, n):
    return H[i, n]

    

def CB1(i, n):
    return 1/2 * ( Q[i+1, n]**2 + Q[i, n]**2 )

def CB2(i, n):
    return teta * Q[i+1, n]

def CB3(i, n):
    return teta * Q[i, n]

def dPdh(i, n):  #replace with pressure change in the node (i,n)
    return 1  #0?


def Ks(i,n):
    return 1 #!!!!!!!!!!!!!


def D(i,n):
    return Ks(i,n) * Rh(i,n)**(2/3) * w[i]

def CC1(i, n):
    return 1/2 * ( D(i+1, n)**2 + D(i, n)**2 )

def CC4(i, n):
    return teta * D(i+1, n) * Ks(i+1, n) * Rh(i+1, n)**(2/3) / 3 * ( 5*B(i+1, n) - 2*Rh(i+1, n)* dPdh(i+1, n) )
                                                                     

def CC5(i, n):
    return CC4(i-1, n)
    return teta * D(i  , n) * Ks(i  , n) * Rh(i  , n)**(2/3) / 3 * ( 5*B(i  , n) - 2*Rh(i  , n)* dPdh(i  , n) )

def CH1(i, n):
    return 1 / 2/g * (1/S(i+1, n) + 1/S(i, n))

def CH4(i, n):
    return -teta/2/g * B(i+1, n) / S(i+1, n)**2

def CH5(i, n):
    return -teta/2/g * B(i  , n) / S(i  , n)**2

def CE2(i, n):
    return 1/2/delt
def CE3(i, n):
    return 1/2/delt
    



def V(i, n):
    return Q[i, n] / ( H[i,n] * w[i] )
    



def CF1(i, n):
    return 1/delx* ( Q[i+1,n]*V(i+1,n) - Q[i,n]*V(i,n) )

def CF2(i, n):
    return 2*teta/delx *V(i+1, n)

def CF3(i, n):
    return -2*teta/delx *V(i, n)

def CH4(i, n):
    return -teta/delx * V(i+1, n)**2 * B(i+1, n)

def CH5(i, n):
    return -teta/delx * V(i  , n)**2 * B(i  , n)
    


def CD1(i,n):
    return 1/2*(h[i+1]-h[i])
def CD4(i,n):
    return teta/delx
def CD5(i,n):
    return -teta/delx

    

def CL(i, n):
    return CB2(i, n) + CC1(i, n) * CH1(i, n) * (CE2(i, n) * CF2(i,n))


def CM(i, n):
    return CC4(i,n) * ( CD1(i,n) + CH1(i,n)*CF1(i,n)) + CC1(i,n)*(CD4(i,n) + CH1(i,n)*CH4(i,n) + CH4(i,n)*CF1(i,n) )

def CN(i, n):
    return -(CB3(i,n) + CC1(i,n)*CH1(i,n) *(CE3(i,n)+CF3(i,n)))

def CO(i, n):
    return -(CC5(i,n)*(CD1(i,n) + CH1(i,n)*CF1(i,n)) + CC1(i,n)*(CD5(i,n)+CH1(i,n)*CH5(i,n)+CH5(i,n)*CF1(i,n)))


def CP(i,n):
    return -(CB1(i,n) + CC1(i,n)* (CD1(i,n) + CH1(i,n)*CF1(i,n)))




delx = 10#0
delt = 1#60

teta=1/2
g=9.81

xmax = 320#00#00
tmax=32*60

x=[]
xx = 0
while xx <= xmax:
    x.append(xx)
    xx += delx

imax = len(x)

t=[]
tt=0
while tt <=tmax:
    t.append(tt)
    tt += delt

#Longitudinal valley base profile , data numeritized from Figure-8
#[longitudinal position, depth]
X_vs_H=[[0, 129], [1350, 117], [3000, 112], [5750, 82], [11000, 60], [13350, 39], [18650, 21], [22300, 18], [32000, 0.1], [50000, 0.1]]

#Positioning and numerization of transversal profiles
# According to section 3.1.4, there are 15 transversal profiles taken in the study from geograpic data.
# However numeritization cannot be done for the cross sections from the information in Figure-2 
# A very rough prediction /approximation may be done from reading the data points in Figure-11
#Width data of Wadi, Figure-11
#[longitudinal position, width of wadi at that position] 
X_vs_W=[[0, 20], [1350, 20], [3000, 100], [5750, 200], [11000, 200], [13350, 2000], [18650, 1500], [22300, 2000], [32000, 3600], [50000, 3000]]

# Upsetream boundry condition (numerization of Figure-6)
# data will be interpolated wrt time discritization scheme  
# [second, m3/s]
t_vs_Q0=[[0,0], [9*60, 15355], [32*60, 0], [1000*60, 0]]



def interpolate_h(x): #interpolates valley height from X_vs_H data
    for i in range(len(X_vs_H)-1):
        xh0 = X_vs_H[i]
        xh1 = X_vs_H[i+1]

        if xh1[0] >= x:
            x0, h0 = xh0[0], xh0[1]
            x1, h1 = xh1[0], xh1[1]
            break
    return h0 + ((x-x0) / (x1-x0)) * (h1-h0) 

h=[]
for i in x:
    h.append(interpolate_h(i))



def interpolate_w(x): #interpolates valley width from X_vs_W data
    for i in range(len(X_vs_W)-1):
        xw0 = X_vs_W[i]
        xw1 = X_vs_W[i+1]

        if xw1[0] >= x:
            x0, w0 = xw0[0], xw0[1]
            x1, w1 = xw1[0], xw1[1]
            break
    return w0 + ((x-x0) / (x1-x0)) * (w1-w0) 

w=[]
for i in x:
    w.append(interpolate_w(i))




def interpolate_Q0(t):
    for i in range(len(t_vs_Q0)-1):
        tq0 = t_vs_Q0[i]
        tq1 = t_vs_Q0[i+1]

        if tq1[0] >= t:
            t0, q0 = tq0[0], tq0[1]
            t1, q1 = tq1[0], tq1[1]
            break
    return q0 + ((t-t0) / (t1-t0))* (q1-q0)


    

H={}
Q={}
for i in range(imax):
    H[i, 0]= h[i]
    Q[i, 0]= 0

#initialization of failure hydogram Figure-6
for n in range(0, tmax//delt):
    Q[0, n] = interpolate_Q0(n*delt)    




for n in range(tmax//delt):
    MAT=[]

    mat = [0]*2*imax
    mat[0] = 1
    MAT.append(mat)


    A=[0]*2*imax
    A[0] = Q[0, n+1] - Q[0, n]   #from failure hydrogram
    A[-1] = 0  #sea level is constant


    mat = [0]*imax
    for i in range(0, imax-1):

        mat = [0]*(2*i) + [-CI(i,n), -CJ(i,n), CG(i,n), CH(i,n)] + [0]*(2*imax-2*i-4)
        A[2*i] = CK(i,n)
        MAT.append(mat)

        mat = [0]*(2*i) + [-CN(i,n), -CO(i,n), CL(i,n), CM(i,n)] + [0]*(2*imax-2*i-4)
        A[2*i+1] = CP(i,n)
        MAT.append(mat)


    mat = [0]*2*imax
    mat[-1] = 1
    MAT.append(mat)


    A2 = matrix_mul2(inverse(MAT), A)

    delQ=[]
    delH=[]
    for j in range(0,len(A2),2):
        delQ.append(A2[j])
        delH.append(A2[j]+1)

    for i in range(imax):
        Q[i, n+1] = Q[i, n] + delQ[i]
        H[i, n+1] = H[i, n] + delH[i]
        

    print(delH[1:20])

