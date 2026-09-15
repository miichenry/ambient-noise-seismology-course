
import numpy, os
import numpy as np
from scipy.signal import fftconvolve, hilbert
from scipy import signal
from DWTO.dwto import DWTO

import matplotlib.pyplot as plt
import psutil

class Correlator:

    fin = []
    wav = []

    cmpname=[]
    ic1=[]
    ic2=[]

    def __init__(self, fin, ccdir, figdir, wav='db1'):
        self.fin = fin
        self.wav = wav
        self.ccdir = ccdir
        self.figdir = figdir

        self.cmpname = ['ZZ', 'ZR', 'ZT', 'RR', 'RT', 'TT']
        self.ic1     = [ 0,    0,    0,    1,    1,    2]
        self.ic2     = [  0,    1,    2,    1,    2,    2]

        if not os.path.isdir(ccdir):
            os.mkdir(ccdir)

        if not os.path.isdir(figdir):
            os.mkdir(figdir)

    def getidxcmp(self, cmp):
        try:
            i=self.cmpname.index(cmp)
        except ValueError:
            return -1, -1
        return self.ic1[i], self.ic2[i]

    def smooth(self, x, sm):
        w = numpy.ones(sm, 'd') / sm
        y = fftconvolve(x, w, mode='same')
        return y

    def whitening(self, x, sm):
        f1 = np.fft.rfft(x)

        if sm!=-1:
            am = np.abs(f1)
            am = self.smooth(am,sm)
            f2 = am * np.exp(1j * np.angle(f1))
        else:
            f2 =np.exp(1j * np.angle(f1))

        return np.fft.irfft(f2)

    def onebit(self, sig, p):
        if p==0:
            return np.sign(sig)
        else:
            return np.abs(sig) ** p * np.sign(sig)

    def stack_linear(self, cc):
        return numpy.sum(cc,0)

    def bandpass(self, fmin, fmax, sig, dt):
        # Filter design
        nyq = 0.5 / dt
        low = fmin / nyq
        high = fmax / nyq
        order = 2
        b, a,  = signal.butter(order, [low, high], btype='band')

        # Windowing
        sig1 = signal.detrend(sig) * signal.tukey(len(sig),0.2)

        # Filtering
        return signal.filtfilt(b, a, sig)

    """
    def stack_pw(self, cc, nu, sm):
        hi = hilbert(cc)
        am = np.sqrt( hi * np.conj(hi) )
        hi=hi/am
        ph=np.sum

        return numpy.sum(cc,0)
    """

    def correlate(self, fmin, fmax, components=['ZZ'], smooth_wh=1, onebit_p=1.):

        import matplotlib
        matplotlib.use('agg')

        test = numpy.zeros(self.fin.nptw)
        testd = DWTO(test, self.wav, dt=self.fin.delta)
        nsc = testd.maxsc

        #log=open('mangosta_cc.log','w')

        sta = self.fin.stations

        for f in range(nsc):

            freq = testd.getfreq(f)
            #print("Wavelet scale ", f, " freq=", freq)
            # log.write("Wavelet scale %d  freq=%f\n" % (f,freq))
            print("Wavelet scale %d  freq=%f" % (f, freq))

            if freq < fmin or freq > fmax:
                print("This frequency is discarded")
                continue

            for cm in range(len(components)):

                for i in range(len(sta)):
                    for j in range(i + 1, len(sta)):
                        nw = self.fin.getnwindows(sta[i], sta[j])

                        #print("Correlating ",sta[i]," with ",sta[j])
                        #log.write("Correlating "+sta[i]+" with "+sta[j]+"\n")
                        #print("Correlating " + sta[i] + " with " + sta[j])

                        print("PAIR ", sta[i], " - ", sta[j]," COMPONENT ",components[cm].upper()," SCALE ",f," FREQ. ",freq)

                        cc = numpy.empty((nw, self.fin.nptw))

                        #print("Created matrix")

                        for t in range(nw):

                            print("Time",self.fin.gettwindow(sta[i], sta[j], t))

                            #log.write("Time %s\n" % self.fin.gettwindow(sta[i], sta[j], t))
                            #print("Time %s" % self.fin.gettwindow(sta[i], sta[j], t))

                            win = self.fin.getwindow(sta[i], sta[j], t)

                            #print("Doing cc on window")
                            ic1, ic2 = self.getidxcmp(components[cm].upper())
                            w1=win[ic1]
                            w2=win[ic2+3]

                            # Spectral whitening
                            if smooth_wh>=3 or smooth_wh==-1:
                                w1 = self.whitening(w1, smooth_wh)
                                w2 = self.whitening(w2, smooth_wh)

                            d1 = DWTO(w1, self.wav, dt=self.fin.delta)
                            d2 = DWTO(w2, self.wav, dt=self.fin.delta)

                            s1 = d1.getsigsc(f)
                            s2 = d2.getsigsc(f)

                            # One-bit normalization
                            if onebit_p<1 and onebit_p>=0:
                                s1 = self.onebit(s1, onebit_p)
                                s2 = self.onebit(s2, onebit_p)

                            cc[t, :] = signal.correlate(s1, s2, mode='same', method='fft')

                            del win, w1, w2, d1, d2, s1, s2

                            print("MEM=", psutil.virtual_memory())

                        filename = "%s/cc_%s_%s_%s_%02d.npy" % (self.ccdir, components[cm].upper(), sta[i], sta[j], f)
                        #print("Writing .npy")
                        numpy.save(filename, cc)

                        #print("Doing stack")
                        st=self.stack_linear(cc)
                        st=st/numpy.std(st)
                        #st=st/numpy.max(numpy.abs(st))

                        #print("del cc")
                        del cc

                        #print("Creating figure")
                        plt.figure()
                        plt.plot(st)
                        filename="%s/figcc_%s_%s_%s_%02d.png" % (self.figdir, components[cm].upper(), sta[i], sta[j], f)
                        plt.savefig(filename)
                        plt.close('all')

                        #print("Write stack")
                        filename="%s/stack_cc_%s_%s_%s_%02d.dat" % (self.ccdir, components[cm].upper(), sta[i], sta[j], f)
                        #print("Output "+filename)
                        numpy.savetxt(filename, st)

                        #print("END LOOP")

    def correlate_trad(self, fmin, fmax, components=['ZZ'], smooth_wh=1, onebit_p=1.):

        sta = self.fin.stations

        for cm in range(len(components)):

            for i in range(len(sta)):
                for j in range(i + 1, len(sta)):
                    nw = self.fin.getnwindows(sta[i], sta[j])

                    print("PAIR ", sta[i], " - ", sta[j]," COMPONENT ",components[cm].upper())

                    cc = numpy.empty((nw, self.fin.nptw))

                    for t in range(nw):

                        print("Time", self.fin.gettwindow(sta[i], sta[j], t))

                        win = self.fin.getwindow(sta[i], sta[j], t)

                        ic1, ic2 = self.getidxcmp(components[cm].upper())
                        w1=win[ic1]
                        w2=win[ic2+3]

                        # Band-pass
                        w1 = self.bandpass(fmin, fmax, w1, self.fin.delta)
                        w2 = self.bandpass(fmin, fmax, w2, self.fin.delta)

                        # Spectral whitening
                        if smooth_wh>=3 or smooth_wh==-1:
                            w1 = self.whitening(w1, smooth_wh)
                            w2 = self.whitening(w2, smooth_wh)

                        # One-bit normalization
                        if onebit_p<1 and onebit_p>=0:
                            w1 = self.onebit(w1, onebit_p)
                            w2 = self.onebit(w2, onebit_p)

                        cc[t, :] = signal.correlate(w1, w2, mode='same', method='fft')

                        del win, w1, w2

                        print("MEM=", psutil.virtual_memory())

                    #filename = "%s/cc_%s_%s_%s_td.npy" % (self.ccdir, components[cm].upper(), sta[i], sta[j])
                    #numpy.save(filename, cc)

                    st=self.stack_linear(cc)
                    st=st/numpy.std(st)

                    del cc

                    plt.figure()
                    plt.plot(st)
                    filename="%s/figcc_%s_%s_%s_td.png" % (self.figdir, components[cm].upper(), sta[i], sta[j])
                    plt.savefig(filename)
                    plt.close('all')

                    filename="%s/stack_cc_%s_%s_%s_td.dat" % (self.ccdir, components[cm].upper(), sta[i], sta[j])
                    numpy.savetxt(filename, st)
