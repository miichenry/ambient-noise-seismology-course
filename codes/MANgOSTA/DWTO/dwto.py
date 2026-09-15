import pywt
from copy import deepcopy
from numpy.fft import rfft, rfftfreq
from numpy import abs, argmax

class DWTO:

    nsam = 0
    wav = 'db1'
    maxsc = 0
    dwta = []

    def __init__(self, sig, wtype='db1', dt=1.0):
        self.wav = wtype
        self.dwta = pywt.wavedec(sig, wtype)
        self.nsc = len(self.dwta)
        self.maxsc = self.nsc-2
        self.nsam = len(sig)

        # Determines empirically the central frequency

        # 1-Set the central coefficient of the 5th scale to one
        s0=5
        n0=self.getncoefsc(s0)
        cd = deepcopy(self.dwta)
        for i in range(self.nsc):
            cd[i] = cd[i] * 0.
            if i == len(self.dwta)-s0-1:
                cd[i][int(n0/2)] = 1.

        sig=pywt.waverec(cd, self.wav)
        #import matplotlib.pyplot as plt
        #plt.plot(sig)
        #plt.show()

        # 2-Compute spectrum and get the central frequency
        sp=abs(rfft(sig))
        fr=rfftfreq(self.nsam, dt)
        imax=argmax(sp)
        self.f0=fr[imax] * 2**s0

        #print("f0=",self.f0)

    def getfreq(self, sc):
        return self.f0 / (2**sc)

    def getcoefsc(self, sc):
        return self.dwta[-sc-1]

    def getncoefsc(self, sc):
        return len(self.dwta[-sc-1])

    def getcoef(self, sc, sh):
        return self.dwta[-sc-1][sh]

    def setcoef(self, sc, sh, val):
        self.dwta[-sc-1][sh] = val

    def getsigsc(self, sc):
        cd=deepcopy(self.dwta)
        for i in range(self.nsc):
            if i!=self.nsc-sc-1:
                cd[i]=cd[i]*0.
        return pywt.waverec(cd, self.wav)
