
import pywt
import numpy as np
from copy import deepcopy

class DWTO2:

    def __init__(self, N, wtype='db1'):

        self.wtype=wtype
        self.base=np.zeros((2**N,2**N))
        self.coef = pywt.wavedec2(self.base, self.wtype)
        self.nscales=len(self.coef)

        self.ncoefsc = []   # Number of coefficients for each scale
        self.nidxsc = []    # Number of axis indices for each scale
        self.idxsc = []     # Size of axis for each scale

        for k in range(self.nscales):
            qq=np.array(self.coef[k])
            self.ncoefsc.append(qq.size)
            self.nidxsc.append(len(qq.shape))
            self.idxsc.append(qq.shape)

        #print(self.ncoefsc)
        #print(self.nidxsc)
        #print(self.idxsc)

    def getnscales(self):
        return self.nscales

    def getncoefscale(self, sc):
        return self.ncoefsc[sc]

    def gettotcoefscale(self, sc):
        return int(np.sum(self.ncoefsc[:sc+1]))

    def getncoef(self):
        return int(np.sum(self.ncoefsc))

    def setcoefsc(self, sc, m, val):

        #if self.nidxsc[sc] == 1:
        #    self.coef[sc][0][0] = val

        if self.nidxsc[sc]==2:
            l = 0
            for i in range(self.idxsc[sc][0]):
                for j in range(self.idxsc[sc][1]):
                    if m == l:
                        self.coef[sc][i][j] = val
                    l = l + 1
        else:
            l = 0
            for i in range(self.idxsc[sc][0]):
                for j in range(self.idxsc[sc][1]):
                    for k in range(self.idxsc[sc][2]):
                        if m == l:
                            self.coef[sc][i][j][k] = val
                        l = l + 1

    def getcoefsc(self, sc, m):

        #print("getcoefsc sc=",sc," m=",m)
        #print("nidxsc=",self.nidxsc[sc])

        #if self.nidxsc[sc] == 1:
        #    return self.coef[sc][0][0]

        if self.nidxsc[sc]==2:
            l = 0
            for i in range(self.idxsc[sc][0]):
                for j in range(self.idxsc[sc][1]):
                    if m == l:
                        return self.coef[sc][i][j]
                    l = l + 1
        else:
            l = 0
            for i in range(self.idxsc[sc][0]):
                for j in range(self.idxsc[sc][1]):
                    for k in range(self.idxsc[sc][2]):
                        if m == l:
                            return self.coef[sc][i][j][k]
                        l = l + 1

    def setcoef(self, m, val):
        mm = m
        sc = 0

        while mm >= self.ncoefsc[sc]:
            mm=mm-self.ncoefsc[sc]
            sc=sc+1

        self.setcoefsc(sc, mm, val)

    def getcoef(self, m):
        mm = m
        sc = 0

        while mm >= self.ncoefsc[sc]:
            mm=mm-self.ncoefsc[sc]
            sc=sc+1

        return self.getcoefsc(sc, mm)

    def getidwt(self):
        self.base = pywt.waverec2(self.coef, self.wtype)
        return self.base

    def sethomoval(self, v0):
        self.base = self.base*0 + v0
        self.coef = pywt.wavedec2(self.base, self.wtype)

        # Reset coefficients for all scales but first
        for sc in range(1,self.nscales):
            nc=self.getncoefscale(sc)
            for j in range(nc):
                self.setcoefsc(sc,j,0.)

    def setval(self, val):
        self.base = val
        self.coef = pywt.wavedec2(self.base, self.wtype)

    def reset(self, value=0.):
        self.base = self.base*0 + value
        self.coef = pywt.wavedec2(self.base, self.wtype)
