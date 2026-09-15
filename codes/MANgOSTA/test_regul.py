
import numpy as np
import matplotlib.pyplot as plt

M=10
xmg, ymg = np.meshgrid( np.linspace(-1,1,M), np.linspace(0,2,M) )
mtrueg=np.zeros((M,M))
mtrueg[3:8,3:7]=1.

xm=np.matrix.flatten(xmg)
ym=np.matrix.flatten(ymg)
mtrue=np.matrix.flatten(mtrueg)

N=50
xd=np.linspace(-2,2,N)
H=0.5
hm2=(ym+H)**2

# Kernel
G=np.zeros((N,M**2))
for i in range(N):
    dist = np.sqrt((xd[i] - xm) ** 2 + hm2 )
    G[i,:] = 1./dist

dtrue=np.dot(G,mtrue) + 0.000 * np.random.randn(N)

# INVERSION

# LSQ
#G1=np.dot(G.T,G)
#G2=np.linalg.inv(G1)
#Gig=np.dot(G2,G.T)

# SVD
#Gig=np.linalg.pinv(G)

# DLSQ
#eps2=10
#G1=np.dot(G.T,G)+np.eye(M**2)*eps2
#G2=np.linalg.inv(G1)
#Gig=np.dot(G2,G.T)

# L-DLSQ
eps2=np.logspace(-10,10,500)
lm=np.empty(500)
ld=np.empty(500)
G0=np.dot(G.T,G)
for i in range(500):
    G1=G0+np.eye(M**2)*eps2[i]
    G2=np.linalg.inv(G1)
    Gig=np.dot(G2,G.T)
    mest = np.dot(Gig, dtrue)
    dsynth = np.dot(G, mest)
    res=dtrue-dsynth
    lm[i]=np.log(np.linalg.norm(mest))
    ld[i]=np.log(np.linalg.norm(res))

curv=np.zeros(500)
for i in range(1,499):
    zeta1 = (ld[i + 1] - ld[i - 1]) / 2
    zeta2 = (ld[i + 1] + ld[i - 1] - 2 * ld[i])
    eta1 = (lm[i + 1] - lm[i - 1]) / 2
    eta2 = (lm[i + 1] + lm[i - 1] - 2 * lm[i])
    curv_num = zeta1 * eta2 - zeta2 * eta1
    curv_den = np.power(zeta1 ** 2 + eta1 ** 2, 1.5)
    curv[i] = curv_num/curv_den

sm=np.array([0.25,0.5,0.75,1,0.75,0.5,0.25])
curv_sm=np.convolve(curv,sm, mode='same')

dlm=np.max(lm)-np.min(lm)
dld=np.max(ld)-np.min(ld)
lm=(lm-np.min(lm))/dlm
ld=(ld-np.min(ld))/dld
dist = np.sqrt(lm**2 + ld**2)
idx = np.argmin(dist)
print("L-DLSQ dist=",idx,eps2[idx])
G1=G0+np.eye(M**2)*eps2[idx]
G2=np.linalg.inv(G1)
Gig=np.dot(G2,G.T)

imax=np.argmax(curv_sm)
print("L-DLSQ curv=",imax,eps2[imax])
G1=G0+np.eye(M**2)*eps2[imax]
G2=np.linalg.inv(G1)
Gig=np.dot(G2,G.T)

m1=np.dot(Gig,dtrue)

# GCV-DLSQ
eps2=np.logspace(-10,10,500)
gcv=np.zeros(500)
G0=np.dot(G.T,G)
for i in range(500):
    G1=G0+np.eye(M**2)*eps2[i]
    G2=np.linalg.inv(G1)
    Gig=np.dot(G2,G.T)
    mest = np.dot(Gig, dtrue)
    dsynth = np.dot(G, mest)
    res=dtrue-dsynth
    A=np.dot(G,Gig)
    tr=np.sum(1.-np.diag(A))**2
    gcv[i]=np.linalg.norm(res)**2/tr

imin=np.argmin(gcv)
print("GCV-DLSQ=",imin,eps2[imin])

"""
plt.plot(gcv,'ko-')
plt.plot(imin,gcv[imin],'ro')
plt.yscale('log')
plt.show()
"""

G1=G0+np.eye(M**2)*eps2[imin]
G2=np.linalg.inv(G1)
Gig=np.dot(G2,G.T)
m2=np.dot(Gig,dtrue)

"""
plt.subplot(1,2,2)
plt.plot(curv_sm,'ko-')
#plt.plot(dist[imin-300:imin+1],'go')
plt.plot(imax,curv_sm[imax],'bo')
#plt.show()
#exit(1)

plt.subplot(1,2,1)
plt.plot(lm,ld,'ko-')
plt.plot(lm[0],ld[0],'ro')
plt.plot(lm[imax],ld[imax],'bo')
plt.show()
#exit(1)
"""

# L-TSVD
U, S, VT = np.linalg.svd(G)
V = VT.T

maxp=len(S)

lm = []
ld = []
lp = []
for p in range(1,maxp):

    Sp = S[range(0, p)]
    Up = U[:, range(0, p)]
    Vp = V[:, range(0, p)]
    iSp = np.diag(1. / Sp)
    Gig = np.dot(np.dot(Vp, iSp), Up.T)
    mest = np.dot(Gig,dtrue)
    sy = np.dot(G,mest)
    res = dtrue-sy

    lm.append(np.log(np.linalg.norm(mest)))
    ld.append(np.log(np.linalg.norm(res)))
    lp.append(p)

# Distance from the origin
dlm=np.max(lm)-np.min(lm)
dld=np.max(ld)-np.min(ld)
lm=(lm-np.min(lm))/dlm
ld=(ld-np.min(ld))/dld
dist = np.sqrt(lm**2 + ld**2)
idx = np.argmin(dist)
popt = lp[idx]
print("L-TSVD=",popt)

"""
plt.subplot(1,2,1)
plt.plot(lm,ld,'ko-')
plt.plot(lm[idx],ld[idx],'ro')

plt.subplot(1,2,2)
plt.plot(dist,'ko-')
plt.plot(idx, dist[idx], 'ro')
plt.show()
"""

Sp = S[range(0, popt)]
Up = U[:, range(0, popt)]
Vp = V[:, range(0, popt)]
iSp = np.diag(1. / Sp)
Gig = np.dot(np.dot(Vp, iSp), Up.T)
m3=np.dot(Gig,dtrue)

# GCV-TSVD
gcv=[]
lp=[]
for p in range(1,maxp):
    Sp = S[range(0, p)]
    Up = U[:, range(0, p)]
    Vp = V[:, range(0, p)]
    iSp = np.diag(1. / Sp)
    Gig = np.dot(np.dot(Vp, iSp), Up.T)
    mest = np.dot(Gig,dtrue)
    sy = np.dot(G,mest)
    res = dtrue-sy
    A = np.dot(G, Gig)
    tr = np.sum(1. - np.diag(A)) ** 2
    gcv.append(np.linalg.norm(res) ** 2 / tr)
    lp.append(p)

imin=np.argmin(gcv)
popt=lp[imin]
print("GCV-TSVD=",popt)
Sp = S[range(0, popt)]
Up = U[:, range(0, popt)]
Vp = V[:, range(0, popt)]
iSp = np.diag(1. / Sp)
Gig = np.dot(np.dot(Vp, iSp), Up.T)
m4=np.dot(Gig,dtrue)

# L-TIKH
L = np.empty((M**2,M**2))
count=0
for i in range(M):
    for j in range(M):
        stencil=np.zeros((M,M))
        stencil[i,j]=-3
        if i > 0: stencil[i - 1, j] = 0.5
        if i < M - 1: stencil[i + 1, j] = 0.5
        if j > 0: stencil[i, j - 1] = 0.5
        if j < M - 1: stencil[i, j + 1] = 0.5
        if i > 0 and j > 0: stencil[i - 1, j - 1] = 0.25
        if i > 0 and j < M - 1: stencil[i - 1, j + 1] = 0.25
        if i < M - 1 and j > 0: stencil[i + 1, j - 1] = 0.25
        if i < M - 1 and j < M - 1: stencil[i + 1, j + 1] = 0.25

        Lrow = np.reshape(stencil, stencil.size)
        idx = np.argwhere(np.abs(Lrow)>0)
        sum = np.sum(Lrow)
        NN = len(idx)
        Lrow[idx] = Lrow[idx] - sum/NN
        L[count,:] = Lrow
        count = count + 1

#plt.pcolor(L)
#plt.show()
eps2=np.logspace(-10,10,500)
lm=np.empty(500)
ld=np.empty(500)
G0=np.dot(G.T,G)
for i in range(500):
    G1=G0+L*eps2[i]
    G2=np.linalg.inv(G1)
    Gig=np.dot(G2,G.T)
    mest = np.dot(Gig, dtrue)
    dsynth = np.dot(G, mest)
    res=dtrue-dsynth
    lm[i]=np.log(np.linalg.norm(mest))
    ld[i]=np.log(np.linalg.norm(res))

curv=np.zeros(500)
for i in range(1,499):
    zeta1 = (ld[i + 1] - ld[i - 1]) / 2
    zeta2 = (ld[i + 1] + ld[i - 1] - 2 * ld[i])
    eta1 = (lm[i + 1] - lm[i - 1]) / 2
    eta2 = (lm[i + 1] + lm[i - 1] - 2 * lm[i])
    curv_num = zeta1 * eta2 - zeta2 * eta1
    curv_den = np.power(zeta1 ** 2 + eta1 ** 2, 1.5)
    curv[i] = curv_num/curv_den

sm=np.array([0.25,0.5,0.75,1,0.75,0.5,0.25])
curv_sm=np.convolve(curv,sm, mode='same')

dlm=np.max(lm)-np.min(lm)
dld=np.max(ld)-np.min(ld)
lm=(lm-np.min(lm))/dlm
ld=(ld-np.min(ld))/dld
dist = np.sqrt(lm**2 + ld**2)
idx = np.argmin(dist)
print("L-TIKH dist=",idx,eps2[idx])
G1=G0+np.eye(M**2)*eps2[idx]
G2=np.linalg.inv(G1)
Gig=np.dot(G2,G.T)

imax=np.argmax(curv_sm)
print("L-TIKH curv=",imax,eps2[imax])
G1=G0+np.eye(M**2)*eps2[imax]
G2=np.linalg.inv(G1)
Gig=np.dot(G2,G.T)

m5=np.dot(Gig,dtrue)

################################33
plt.subplot(2,3,1)
plt.pcolor(xmg,ymg,np.resize(m1,(M,M)))
plt.gca().invert_yaxis()
plt.gca().title.set_text('L-DLSQ')

plt.subplot(2,3,4)
plt.pcolor(xmg,ymg,np.resize(m2,(M,M)))
plt.gca().invert_yaxis()
plt.gca().title.set_text('GCV-DLSQ')

plt.subplot(2,3,2)
plt.pcolor(xmg,ymg,np.resize(m3,(M,M)))
plt.gca().invert_yaxis()
plt.gca().title.set_text('L-TSVD')

plt.subplot(2,3,5)
plt.pcolor(xmg,ymg,np.resize(m4,(M,M)))
plt.gca().invert_yaxis()
plt.gca().title.set_text('GCV-TSVD')

plt.subplot(2,3,3)
plt.pcolor(xmg,ymg,np.resize(m4,(M,M)))
plt.gca().invert_yaxis()
plt.gca().title.set_text('L-TIKH')

plt.show()

exit(1)

############################
# OUTPUT

mest=np.dot(Gig,dtrue)
dsynth=np.dot(G,mest)

plt.subplot(2,2,1)
plt.plot(xd,dtrue,'rx-')
plt.subplot(2,2,3)
plt.pcolor(xmg,ymg,mtrueg)
plt.gca().invert_yaxis()
plt.xlim([-2,2])
plt.ylim([-1,2])

plt.subplot(2,2,2)
plt.plot(xd,dtrue,'rx-')
plt.plot(xd,dsynth,'kx-')
plt.subplot(2,2,4)
plt.pcolor(xmg,ymg,np.resize(mest,(M,M)))
plt.gca().invert_yaxis()
plt.xlim([-2,2])
plt.ylim([-1,2])

plt.show()
