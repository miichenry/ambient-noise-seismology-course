
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider, TextBox
from copy import deepcopy
import sys, os
import tkinter

class FTAN:

    def __init__(self):

        # Init parameters
        self.P1 = 1
        self.P2 = 10
        self.NPER = 100
        self.PERIODS = np.logspace(np.log10(self.P1), np.log10(self.P2), self.NPER)

        self.V1 = 0.5
        self.V2 = 4
        self.NVEL = 100
        self.VEL = np.linspace(self.V1, self.V2, self.NVEL)

        self.curp = []
        self.curv = []

        self.newp = []
        self.newv = []

        self.flagdisp = False

        self.curfile = ''

        ##################################3
        """
        # Read data
        self.sig = np.loadtxt('stack_cc_ZZ_TCHA_TVIL_td.dat')
        fft = np.fft.rfft(self.sig)
        wfft = np.exp(1j * np.angle(fft))
        self.sig = np.fft.irfft(wfft)
        self.sig = self.butter_bandpass(self.sig, 1. / self.P2, 1. / self.P1, 100.)
        self.dist = 9
        #self.sig = self.sig * np.hanning(len(self.sig))
        """

        ##################################3
        """
        # Synth data
        NPT = 4000
        DELTA = 0.01
        T = np.arange(0, NPT) * DELTA - (NPT*DELTA)/2
        sig = T * 0
        self.dist = 15

        # Synthetic
        for i in range(1000):
            lper = np.random.uniform(-1, 0)
            per = 10 ** lper
            freq = 1 / per
            vel = 2. * (lper + 1) + 1
            # print lper, per, freq, vel
            t0 = self.dist / vel
            wi = per * 2
            sig = sig + np.exp(-((T - t0) / wi) ** 2) * np.cos(2 * np.pi * freq * (T - t0))
            sig = sig + np.exp(-((T + t0) / wi) ** 2) * np.cos(2 * np.pi * freq * (T + t0))

        # SW
        fft = np.fft.rfft(sig)
        wfft = np.exp(1j * np.angle(fft))
        self.sig = np.fft.irfft(wfft)
        # Taper
        self.sig = self.sig*np.hanning(NPT)
        """

        ##################################3
        NPT1=1000
        DELTA = 0.01
        T = np.arange(0, NPT1) * DELTA
        s1=scipy.signal.chirp(T,1,10,10, phi=90, method='quadratic')
        s1b=np.hstack((np.zeros(500), s1, np.zeros(500)))

        s2b=np.flip(s1b, axis=0)
        self.sig=np.hstack((s2b,s1b))
        self.dist = 15

        ##################################3

        self.sig_ff=np.empty(0)

        self.NPT = len(self.sig)
        self.DELTA = 0.01

        TMAX = self.NPT * self.DELTA / 2
        self.TSISMO = np.arange(0, self.NPT) * self.DELTA - TMAX

        # Init fig
        self.fig, ax = plt.subplots()

        # FTAN
        plt.subplots_adjust(left=0.25, bottom=0.5)
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
        ax_p1 = plt.axes([0.30, 0.43, 0.25, 0.05], facecolor='g')
        self.sl_p1 = Slider(ax_p1, 'log P1', -2, 2, valinit=0)
        ax_p2 = plt.axes([0.63, 0.43, 0.25, 0.05], facecolor='g')
        self.sl_p2 = Slider(ax_p2, 'log P2', -2, 2, valinit=1)

        ax_v1 = plt.axes([0.30, 0.37, 0.25, 0.05], facecolor='g')
        self.sl_v1 = Slider(ax_v1, 'V1', 0.1, 5, valinit=0.5)
        ax_v2 = plt.axes([0.63, 0.37, 0.25, 0.05], facecolor='g')
        self.sl_v2 = Slider(ax_v2, 'V2', 0.1, 5, valinit=4)

        ax_tmax = plt.axes([0.30, 0.31, 0.58, 0.05], facecolor='g')
        self.sl_tmax = Slider(ax_tmax, 'TMAX', 0, 100, valinit=100)

        # SIG
        self.sig_ax = plt.axes([0.25, 0.1, 0.65, 0.19])
        self.psig, = plt.plot(self.TSISMO, self.sig, lw=0.2)
        plt.gca().set_xlabel('Time (s)')
        plt.gca().set_yticks([])
        self.ax_sig=plt.gca()
        plt.xlim([np.min(self.TSISMO), np.max(self.TSISMO)])

        # WIDGETS
        c1_ax = plt.axes([0.05, 0.85, 0.15, 0.05])
        self.c1_tb = TextBox(c1_ax, 'c\u2081 ', '0.0')

        w_ax = plt.axes([0.05, 0.80, 0.15, 0.05])
        self.w_tb = TextBox(w_ax, 'W ', '5.0')

        alpha_ax = plt.axes([0.05, 0.75, 0.15, 0.05])
        self.alpha_tb = TextBox(alpha_ax, '\u03B1 ', '50.0')

        #
        cc_ax = plt.axes([0.05, 0.65, 0.15, 0.05])
        self.cc_tb = TextBox(cc_ax, 'CC', 'cc')

        dd_ax = plt.axes([0.05, 0.60, 0.15, 0.05])
        self.dd_tb = TextBox(dd_ax, 'DISP', 'disp')

        dfile_ax = plt.axes([0.05, 0.55, 0.15, 0.05])
        self.dfile_tb = TextBox(dfile_ax, 'DFILE', 'pair_distance.dat')

        #
        cur_ax = plt.axes([0.05, 0.45, 0.15, 0.05])
        self.cur_tb = TextBox(cur_ax, 'FILE', '')

        dist_ax = plt.axes([0.05, 0.40, 0.15, 0.05])
        self.dist_tb = TextBox(dist_ax, 'DIST', '')

        ref_ax = plt.axes([0.05, 0.35, 0.15, 0.05])
        self.ref_tb = TextBox(ref_ax, 'REF', 'disp_RR.dat')

        #
        ap_ax = plt.axes([0.05, 0.25, 0.15, 0.05])
        ap_bt = Button(ap_ax, 'Apply')
        ap_bt.on_clicked(self.update_all)

        res_ax = plt.axes([0.05, 0.20, 0.15, 0.05])
        res_bt = Button(res_ax, 'Reset')
        res_bt.on_clicked(self.reset)

        op_ax = plt.axes([0.05, 0.15, 0.15, 0.05])
        op_bt = Button(op_ax, 'Open')
        op_bt.on_clicked(self.open)

        sa_ax = plt.axes([0.05, 0.10, 0.15, 0.05])
        sa_bt = Button(sa_ax, 'Save')
        sa_bt.on_clicked(self.save)

        self.update_all(None)

        plt.show()

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

        self.P1 = 10**self.sl_p1.val
        self.P2 = 10**self.sl_p2.val
        self.NPER = 100
        self.PERIODS = np.logspace(np.log10(self.P1), np.log10(self.P2), self.NPER)

        self.V1 = self.sl_v1.val
        self.V2 = self.sl_v2.val
        self.NVEL = 100
        self.VEL = np.linspace(self.V1, self.V2, self.NVEL)

        #print(self.P1, self.P2, self.V1, self.V2)

        c1 = float(self.c1_tb.text)
        w = float(self.w_tb.text)
        alpha = float(self.alpha_tb.text)

        self.curp = deepcopy(self.newp)
        self.curv = deepcopy(self.newv)

        if len(self.curp)>1:
            self.curp, self.curv = self.adjcur(self.curp, self.curv)

        self.newp = []
        self.newv = []

        # Spectral whitening, filtering and windowing
        fft = np.fft.rfft(self.sig)
        wfft = np.exp(1j * np.angle(fft))
        sig2 = np.fft.irfft(wfft)
        sig2 = self.butter_bandpass(sig2, 1. / self.P2, 1. / self.P1, 100.)
        sig2 = sig2 * np.hanning(len(sig2))

        freq = np.fft.rfftfreq(self.NPT, self.DELTA)
        ft = np.fft.rfft(sig2)

        # Floating filter
        #print("FF")
        if len(self.curp) > 1:
            pp = np.array(self.curp)
            vv = np.array(self.curv)

            p0 = np.min(pp)
            p1 = np.max(pp)

            f0 = 1./p1
            f1 = 1./p0

            pick_freq = 2 * np.pi / pp
            w0 = 2 * np.pi * f0
            w1 = 2 * np.pi * f1
            int_freq = np.linspace(w0, w1, 100)
            dw = (w1 - w0) / 100
            int_u = np.interp(int_freq, pick_freq, vv)
            integr = np.sum(1 / int_u) * dw
            #print(integr)
            psi = self.dist * integr + c1 * freq * 2 * np.pi

            fact = np.exp(np.complex(0, 1) * psi)
            ft_ff = ft * fact

            sig2 = np.fft.irfft(ft_ff)

        # FTAN
        #print("FTAN")
        AM = np.zeros((self.NPER, self.NVEL))
        ft = np.fft.rfft(sig2)

        #plt.figure()

        for i in range(self.NPER):

            pe = self.PERIODS[i]
            fr = 1. / pe
            #alpha = max(61.6 - 0.8 * pe, 5)
            fgauss = np.exp(- alpha * ((freq - fr) / fr) ** 2)

            ft_filt = ft * fgauss
            sig_filt = np.fft.irfft(ft_filt)

            data_envelope = np.abs(scipy.signal.hilbert(sig_filt))
            #plt.plot(self.TSISMO, data_envelope)

            # Gaussian tapering

            if len(self.curp) > 1:
                u = np.interp(pe, pp, vv)
                tg = self.dist / u
                mask = np.exp( - ((self.TSISMO - tg) / w) ** 2) + np.exp( - ((self.TSISMO + tg) / w) ** 2)
                data_envelope = data_envelope * mask

            for j in range(self.NVEL):
                tg = self.dist / self.VEL[j]
                #print(j, self.VEL[j], tg)
                AM[i, j] = np.interp(tg, self.TSISMO, data_envelope) + np.interp(-tg, self.TSISMO, data_envelope)
                #AM[i, j] = np.interp(tg, self.TSISMO, data_envelope)

            sm = np.max(AM[i, :])
            AM[i, :] = AM[i, :] / sm

        #plt.show()

        # Plot FTAN
        #print("Plot FTAN")
        for coll in self.cftan.collections:
            self.axp.collections.remove(coll)
        self.cftan = self.axp.contourf(self.PERIODS, self.VEL, AM.T, 20)

        #print("Plot MAX")
        imax=np.argmax(AM, axis=1)
        try:
            self.pl_max.set_data(self.PERIODS, self.VEL[imax])
            #print("MAX ydata")
        except AttributeError:
            self.pl_max, = self.axp.plot(self.PERIODS, self.VEL[imax], 'kx')
            #print("MAX plot")

        self.ax_ftan.set_xlim([self.P1, self.P2])
        self.ax_ftan.set_ylim([self.V1, self.V2])

        # Read and plot REFERENCE disp curve
        refdc_file = self.ref_tb.text

        if os.path.isfile(refdc_file):
            refdc = np.loadtxt(refdc_file)
            try:
                self.refdisp1.set_data(refdc[:,0], refdc[:,1])
                self.refdisp2.set_data(refdc[:,0], refdc[:,1])

            except AttributeError:
                self.refdisp1, = self.axp.plot(refdc[:,0], refdc[:,1], 'bo')
                self.refdisp2, = self.axp.plot(refdc[:,0], refdc[:,1], 'b')

        #print("Plot NEWP")
        try:
            self.rdisp1.set_data([],[])
            self.rdisp2.set_data([],[])
            #print("NEWP data")
        except:
            #print("NEWP no")
            pass

        #print("Plot CURP")
        if len(self.curp)>0:

            try:
                self.pdisp1.set_data(self.curp, self.curv)
                self.pdisp2.set_data(self.curp, self.curv)
                #print("CURP data")

            except AttributeError:
                self.pdisp1, = self.axp.plot(self.curp, self.curv, 'ko')
                self.pdisp2, = self.axp.plot(self.curp, self.curv, 'k')
                #print("CURP plot")

            except:
                print(sys.exc_info()[0])

        # Plot SISMO
        #print("Plot SISMO")
        self.psig.set_data(self.TSISMO, sig2)
        mm=np.max(sig2)-np.min(sig2)
        self.sig_ax.set_ylim([np.min(sig2)-mm*0.1, np.max(sig2)+mm*0.1])

        tmaxl=self.sl_tmax.val*np.max(self.TSISMO)/100.
        self.ax_sig.set_xlim([-tmaxl, tmaxl])

        #print("REDRAW")
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

    def getdist(self, name, pairfile):

        P1=name[3:7]
        P2=name[8:12]
        PAIR='%s_%s' % (P1, P2)

        if not os.path.exists(pairfile):
            print('Error pair file %s not found!' % pairfile)
            exit(1)

        num_lines = sum(1 for line in open(pairfile))

        ifile=open(pairfile,'r')
        for i in range(num_lines):
            row=ifile.readline().split()
            if row[0]==PAIR:
                dist=float(row[1])/1000.
                print(dist)
                return dist

        print('Error pair %s not found in file %s!' % (PAIR, pairfile))
        exit(1)

    def load(self, tkwin, listbox):

        sel=listbox.curselection()
        if len(sel)>0:

            #print(self.filelist[sel[0]])
            self.curname=self.namelist[sel[0]]
            self.sig = np.loadtxt(self.filelist[sel[0]])

            """
            fft = np.fft.rfft(self.sig)
            wfft = np.exp(1j * np.angle(fft))
            self.sig = np.fft.irfft(wfft)
            self.sig = self.butter_bandpass(self.sig, 1. / self.P2, 1. / self.P1, 100.)
            """

            self.cur_tb.set_val(self.curname)

            self.dist = self.getdist(self.curname, self.dfile_tb.text)
            self.dist_tb.set_val(self.dist)

            self.NPT = len(self.sig)
            self.DELTA = 0.01

            TMAX = self.NPT * self.DELTA / 2
            self.TSISMO = np.arange(0, self.NPT) * self.DELTA - TMAX
            #print(self.NPT, TMAX)

            tkwin.destroy()
            #del tkwin
            self.update_all(None)

        else:
            tkwin.destroy()
            #del tkwin
            #print("No selection")

    def open(self, val):
        master = tkinter.Tk()

        listbox = tkinter.Listbox(master)
        listbox.pack()

        # check files
        self.filelist=[]
        self.namelist=[]

        ccdir=self.cc_tb.text
        dddir=self.dd_tb.text

        for root, dirnames, filenames in os.walk(ccdir):
            for filename in filenames:
                name=filename[9:24]
                self.filelist.append('%s/%s' % (root,filename))
                self.namelist.append(name)
                dispname='%s/disp_%s.dat' % (dddir, name)
                #print(dispname)
                if os.path.isfile(dispname):
                    listname='(*) %s' % name
                else:
                    listname = '( ) %s' % name
                listbox.insert(tkinter.END, listname)

        b=tkinter.Button(master, text="OK", command=lambda: self.load(master, listbox))
        b.pack()

        tkinter.mainloop()

    def save(self, val):

        dddir = self.dd_tb.text

        # Check dispdir
        if not os.path.exists(dddir): os.makedirs(dddir)

        # Write file
        dispname = '%s/disp_%s.dat' % (dddir, self.curname)
        ofile=open(dispname, 'w')
        for i in range(len(self.curp)):
            ofile.write('%f %f -1\n' % (self.curp[i], self.curv[i]))
        ofile.close()

#############################################
if __name__ == '__main__':

    f = FTAN()
