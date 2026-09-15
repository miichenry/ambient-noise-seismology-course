
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons, Slider, TextBox
from copy import deepcopy
import sys, os

class MFTAN:

    def __init__(self):

        # Init parameters
        self.P1 = 1
        self.P2 = 10
        self.NPER = 10
        self.PERIODS = np.logspace(np.log10(self.P1), np.log10(self.P2), self.NPER)

        self.V1 = 0.5
        self.V2 = 4
        self.NVEL = 10
        self.VEL = np.linspace(self.V1, self.V2, self.NVEL)

        self.curp = []
        self.curv = []

        self.newp = []
        self.newv = []

        self.flagdisp = False

        self.curfile = ''

        ##################################3

        # Init fig
        self.fig, ax = plt.subplots()

        # FTAN
        plt.subplots_adjust(left=0.25, bottom=0.20)
        self.axp = plt.gca()
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        mat = np.zeros((self.NPER, self.NVEL))
        self.cftan = plt.contourf(self.PERIODS, self.VEL, mat.T, 20)
        plt.gca().xaxis.tick_top()
        plt.gca().yaxis.tick_right()
        plt.gca().set_xlabel('Period (s)')
        plt.gca().xaxis.set_label_position("top")
        plt.gca().set_ylabel('Group velocity (km/s)')
        plt.gca().yaxis.set_label_position("right")
        plt.gca().set_xscale('log')
        plt.grid(True)
        self.ax_ftan=plt.gca()

        # SLIDERS
        ax_p1 = plt.axes([0.30, 0.13, 0.25, 0.05], facecolor='g')
        self.sl_p1 = Slider(ax_p1, 'log P1', -2, 2, valinit=0)
        ax_p2 = plt.axes([0.63, 0.13, 0.25, 0.05], facecolor='g')
        self.sl_p2 = Slider(ax_p2, 'log P2', -2, 2, valinit=1)

        ax_v1 = plt.axes([0.30, 0.05, 0.25, 0.05], facecolor='g')
        self.sl_v1 = Slider(ax_v1, 'V1', 0.1, 5, valinit=0.5)
        ax_v2 = plt.axes([0.63, 0.05, 0.25, 0.05], facecolor='g')
        self.sl_v2 = Slider(ax_v2, 'V2', 0.1, 5, valinit=4)

        # TEXT BOX
        alpha_ax = plt.axes([0.05, 0.85, 0.15, 0.05])
        self.alpha_tb = TextBox(alpha_ax, '\u03B1 ', '50.0')

        # RADIO_BUTTONS
        ax_rb = plt.axes([0.05, 0.50, 0.15, 0.30])
        self.rb = RadioButtons(ax_rb, ('RR', 'ZZ', 'TT'))

        # TEXT BOX
        cc_ax = plt.axes([0.05, 0.40, 0.15, 0.05])
        self.cc_tb = TextBox(cc_ax, 'CC', 'cc')

        dd_ax = plt.axes([0.05, 0.35, 0.15, 0.05])
        self.dd_tb = TextBox(dd_ax, 'DISP', 'disp')

        dfile_ax = plt.axes([0.05, 0.30, 0.15, 0.05])
        self.dfile_tb = TextBox(dfile_ax, 'DFILE', 'pair_distance.dat')

        # BUTTON
        ap_ax = plt.axes([0.05, 0.20, 0.15, 0.05])
        ap_bt = Button(ap_ax, 'Apply')
        ap_bt.on_clicked(self.update_all)

        res_ax = plt.axes([0.05, 0.15, 0.15, 0.05])
        res_bt = Button(res_ax, 'Reset')
        res_bt.on_clicked(self.reset)

        sa_ax = plt.axes([0.05, 0.05, 0.15, 0.05])
        sa_bt = Button(sa_ax, 'Save')
        sa_bt.on_clicked(self.save)

        self.update_all(None)

        plt.show()

    def readdata(self):

        print("READ DATA")

        # read dist
        pairfile = self.dfile_tb.text
        num_lines = sum(1 for line in open(pairfile))
        f=open(pairfile,'r')
        pair=[]
        dist=[]
        for i in range(num_lines):
            row=f.readline().split()
            pair.append(row[0])
            dist.append(float(row[1]))
        f.close()

        comp = self.rb.value_selected

        dir = self.cc_tb.text
        minlen=sys.maxsize
        datalist=[]
        distlist=[]
        for root, dirnames, filenames in os.walk(dir):
            for filename in filenames:
                #print(root,filename)
                compfile=filename[9:11]
                if comp==compfile:
                    pairfile=filename[12:21]
                    try:
                        idx=pair.index(pairfile)
                        distlist.append(dist[idx])
                        data=np.loadtxt('%s/%s' % (root,filename))
                        if len(data)<minlen: minlen=len(data)
                        datalist.append(data)
                    except: pass

        nfiles=len(distlist)
        self.datamat=np.empty((nfiles,minlen))
        for i in range(nfiles):
            self.datamat[i,:] = datalist[i][0:minlen]

        # Remove mean
        for i in range(nfiles):
            self.datamat[i,:] = self.datamat[i,:] - np.mean(self.datamat[i,:])

        self.dist = np.array(distlist)
        self.NPT = minlen
        self.DELTA = 0.01

        TMAX = self.NPT * self.DELTA / 2
        self.TSISMO = np.arange(0, self.NPT) * self.DELTA - TMAX

    def calcFTAN(self):

        print("FTAN")

        self.P1 = 10**self.sl_p1.val
        self.P2 = 10**self.sl_p2.val
        self.NPER = 100
        self.PERIODS = np.logspace(np.log10(self.P1), np.log10(self.P2), self.NPER)

        self.V1 = self.sl_v1.val
        self.V2 = self.sl_v2.val
        self.NVEL = 100
        self.VEL = np.linspace(self.V1, self.V2, self.NVEL)

        self.AM = np.zeros((self.NPER, self.NVEL))

        alpha = float(self.alpha_tb.text)
        freq = np.fft.rfftfreq(self.NPT, self.DELTA)

        # Spectral whitening, filtering and windowing
        data=deepcopy(self.datamat)
        for i in range(data.shape[0]):
            fft = np.fft.rfft(data[i,:])
            wfft = np.exp(1j * np.angle(fft))
            sig2 = np.fft.irfft(wfft)
            sig2 = self.butter_bandpass(sig2, 1. / self.P2, 1. / self.P1, 100.)
            data[i,:] = sig2 * np.hanning(len(sig2))

        """
        if len(self.curp) > 1:
            pp = np.array(self.curp)
            vv = np.array(self.curv)
        """

        for i in range(self.NPER):

            pe = self.PERIODS[i]
            fr = 1. / pe
            fgauss = np.exp(- alpha * ((freq - fr) / fr) ** 2)

            # Compute envelopes for the band
            envel = data * 0
            for k in range(data.shape[0]):
                ft = np.fft.rfft(data[k,:])
                ft_filt = ft * fgauss
                sig_filt = np.fft.irfft(ft_filt)
                envel[k,:] = np.abs(scipy.signal.hilbert(sig_filt))

            # FTAN
            for j in range(self.NVEL):
                for k in range(data.shape[0]):
                    tg = self.dist[k] / (self.VEL[j]*1000.)
                    self.AM[i, j] = self.AM[i, j] + np.interp(tg, self.TSISMO, envel[k,:]) + np.interp(-tg, self.TSISMO, envel[k,:])

            # Normalize
            sm = np.max(self.AM[i, :])
            self.AM[i, :] = self.AM[i, :] / sm

    def adjcur(self, curp, curv):

        i=1
        while i<len(curp) and curp[i]>curp[i-1]:
            i=i+1

        if i<len(curp):
            curp = curp[0:i]
            curv = curv[0:i]

        return curp, curv

    def onclick(self, event):

        per = event.xdata
        vel = event.ydata
        #print(per,vel)

        if per!=None and vel!=None and per>=self.P1 and per<=self.P2 and vel>=self.V1 and vel<=self.V2:
            self.newp.append(per)
            self.newv.append(vel)
            #print("CLICK ",per,vel)

            try:
                self.rdisp1.set_data(self.newp, self.newv)
                self.rdisp2.set_data(self.newp, self.newv)

            except AttributeError:
                self.rdisp1, = self.axp.plot(self.curp, self.curv, 'ro')
                self.rdisp2, = self.axp.plot(self.curp, self.curv, 'r')

            self.fig.canvas.draw_idle()

    def butter_bandpass(self, data, lowcut, highcut, fs, order=2):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        if high>0.9: high=0.9
        b, a = scipy.signal.butter(order, [low, high], btype='band')
        y = scipy.signal.lfilter(b, a, data)

        return y

    def update_all(self, val):

        self.readdata()
        self.calcFTAN()

        self.curp = deepcopy(self.newp)
        self.curv = deepcopy(self.newv)

        if len(self.curp)>1:
            self.curp, self.curv = self.adjcur(self.curp, self.curv)

        self.newp = []
        self.newv = []

        # Plot FTAN
        #print("Plot FTAN")
        for coll in self.cftan.collections:
            self.axp.collections.remove(coll)
        self.cftan = self.axp.contourf(self.PERIODS, self.VEL, self.AM.T, 20)

        #print("Plot MAX")
        imax=np.argmax(self.AM, axis=1)
        try:
            self.pl_max.set_data(self.PERIODS, self.VEL[imax])
            #print("MAX ydata")
        except AttributeError:
            self.pl_max, = self.axp.plot(self.PERIODS, self.VEL[imax], 'kx')
            #print("MAX plot")

        self.ax_ftan.set_xlim([self.P1, self.P2])
        self.ax_ftan.set_ylim([self.V1, self.V2])

        #print("Plot NEWP")
        try:
            self.rdisp1.set_data([],[])
            self.rdisp2.set_data([],[])
            #print("NEWP data")
        except:
            pass

        if len(self.curp)>0:

            try:
                self.pdisp1.set_data(self.curp, self.curv)
                self.pdisp2.set_data(self.curp, self.curv)

            except AttributeError:
                self.pdisp1, = self.axp.plot(self.curp, self.curv, 'ko')
                self.pdisp2, = self.axp.plot(self.curp, self.curv, 'k')

            except:
                print(sys.exc_info()[0])

        self.fig.canvas.draw_idle()

    def reset(self, val):

        try:
            self.pdisp1.set_data([],[])
            self.pdisp2.set_data([],[])
        except:
            pass

        self.curp=[]
        self.curv=[]
        self.newp=[]
        self.newv=[]
        self.update_all(None)

    def save(self, val):

        comp = self.rb.value_selected

        # Write file
        dispname = 'disp_%s.dat' % comp
        ofile=open(dispname, 'w')
        for i in range(len(self.curp)):
            ofile.write('%f %f -1\n' % (self.curp[i], self.curv[i]))
        ofile.close()

#############################################
if __name__ == '__main__':

    f = MFTAN()
