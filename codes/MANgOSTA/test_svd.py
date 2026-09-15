
import numpy as np
import matplotlib.pyplot as plt

def pprint(a):
    for n in a: print("%.3f " % n, end='')
    print('')

N=20
M=5

G=np.random.randn(N,M)
#m=np.random.randn(M)
m=np.zeros(M)
m[2]=1.
d=np.dot(G,m)

# PINV
Gig=np.linalg.pinv(G)
mest1=np.dot(Gig,d)

# SVD mat
U, S, VT = np.linalg.svd(G)
V = VT.T
UT = U.T

p=M
Sp = S[range(0, p)]
Up = U[:, range(0, p)]
Vp = V[:, range(0, p)]
iSp = np.diag(1. / Sp)
Gig = np.dot(np.dot(Vp, iSp), Up.T)
mest2=np.dot(Gig,d)

# SVD comp
mest3=m*0
for i in range(M):
    Ui = U[:,i]
    Vi = V[:,i]
    fac = 1
    mest3 = mest3 + fac*np.dot(Ui,d)*Vi/Sp[i]

# DLSQ
alpha = 1.0e+0
G1 = np.dot(G.T,G) + alpha**2 * np.eye(M)
G2 = np.linalg.inv(G1)
Gig = np.dot(G2,G.T)
mest4=np.dot(Gig,d)

# SVD impl of DLSQ
mest5=m*0
for i in range(M):
    Ui = U[:,i]
    Vi = V[:,i]
    fac = Sp[i]**2 / (Sp[i]**2 + alpha**2)
    mest5 = mest5 + fac*np.dot(Ui,d)*Vi/Sp[i]

# Tikhonov
L=np.zeros((M,M))
L[0,0:2]=np.array([-2,1])
#L[0,0:3]=np.array([1,-2,1])
for i in range(1,M-1): L[i,i-1:i+2]=np.array([1,-2,1])
L[M-1,M-2:M]=np.array([1,-2])
#L[M-1,M-3:M]=np.array([1,-2,1])
L=-L

"""
print(L)

count=0
for beta in np.logspace(-1,1,10):
    G1 = np.dot(G.T,G) + beta**2 * L
    G2 = np.linalg.inv(G1)
    Gig = np.dot(G2,G.T)
    mest6=np.dot(Gig,d)
    mest6=mest6/np.max(np.abs(mest6))
    plt.plot(mest6+count,'k-')
    count=count+1
plt.show()
exit(1)
"""

# Tikhonov with GSVD
# mest 7

# Levenberg-Marquandt
gamma=1.0e-0
G1 = np.dot(G.T,G)
DLM = np.diag(G1)
G2 = G1 + gamma**2 * np.diag(DLM)
G3 = np.linalg.inv(G2)
Gig = np.dot(G3,G.T)
mest8=np.dot(Gig,d)

############################
pprint(m)
#pprint(mest1)
#pprint(mest2)
#pprint(mest3)
#pprint(mest4)
#pprint(mest5)
#pprint(mest6)
##pprint(mest7)
pprint(mest8)
