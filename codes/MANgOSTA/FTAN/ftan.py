
import glob, os
import numpy
import obspy
import matplotlib.pyplot as plt

class FTAN:

    ccdir = []
    st = []

    def __init__(self, ccdir, figdir, st):
        self.ccdir = ccdir
        self.figdir = figdir
        self.st = st

    def whitening(self, sig):
        ft=numpy.fft.rfft(sig)
        amp=numpy.abs(ft)
        for i in range(len(amp)):
            if amp[i]>0:
                ft[i]=ft[i]/amp[i]
            else:
                ft[i]=0.0
        return numpy.fft.irfft(ft)

    def plotftan(self, sig, dist, ofile):

        PERIODS = numpy.logspace(numpy.log10(5), numpy.log10(250), 500)
        alpha = 0.02
        NPER = len(PERIODS)

        VEL = numpy.linspace(1, 4, 500)
        NVEL = len(VEL)

        trz=obspy.Trace(sig)
        trz.stats.delta=0.01

        NPT = sig.size
        #DELTA = trz.stats.delta
        DELTA=0.01
        TMAX = NPT * DELTA / 2
        TSISMO = numpy.arange(0, NPT)*DELTA - TMAX

        AM = numpy.zeros((NPER, NVEL))

        for i in range(NPER):

            pe = PERIODS[i]
            fr = 1. / pe

            f1 = fr - alpha * fr
            f2 = fr + alpha * fr

            tr = trz.copy()
            tr.filter('bandpass', freqmin=f1, freqmax=f2, corners=4)
            data_envelope = obspy.signal.filter.envelope(tr.data)

            for j in range(NVEL):
                tg = dist / VEL[j]
                AM[i, j] = numpy.interp(tg, TSISMO, data_envelope) + numpy.interp(-tg, TSISMO, data_envelope)

            # sm=np.sum(AM[i,:])
            sm = numpy.max(AM[i, :])
            AM[i, :] = AM[i, :] / sm

        # AM=np.log(AM.T)
        plt.matshow(AM.T)
        plt.savefig(ofile)
        plt.close('all')
        #plt.show()

    def doFTAN(self, fig=False):

        pairs=[]
        ff=[]
        cmp=[]

        olddir=os.getcwd()
        os.chdir(self.ccdir)
        for file in glob.glob("stack*"):
            ff.append(file)
            cmp.append(file[9:11])
            pairs.append(file[12:21])

        upairs=list(set(pairs))
        ucmp=list(set(cmp))

        for p in range(len(upairs)):
            for c in range(len(ucmp)):

                print("FTAN for comp ",ucmp[c]," station pair ",upairs[p])

                idx=[]
                for i in range(len(ff)):
                    if pairs[i]==upairs[p] and cmp[i]==ucmp[c]:
                        idx.append(i)

                s0=numpy.empty(0)
                for i in range(len(idx)):
                    s1=numpy.loadtxt(ff[idx[i]])
                    s1=s1/numpy.std(s1)

                    #plt.plot(s1)
                    #plt.title('s1')
                    #plt.show()

                    if s0.size==0: s0=s1
                    else: s0=s0+s1

                #plt.plot(s0)
                #plt.title('s0 before')
                #plt.show()

                #s0=self.whitening(s0)
                s0=s0/numpy.std(s0)

                #plt.plot(s0)
                #plt.title('s0 after')
                #plt.show()

                sta1 = upairs[p][0:4]
                sta2 = upairs[p][5:9]
                dist=self.st.getdist(sta1, sta2)

                ofile="%s/%s/fig_ftan_%s_%s.png" % (olddir, self.figdir, ucmp[c], upairs[p])
                print(ofile)
                self.plotftan(s0,dist,ofile)
