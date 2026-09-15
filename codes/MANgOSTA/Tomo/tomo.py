
import numpy as np
import os, fnmatch, utm

from copy import deepcopy
from Tomo.velmod import Velmod
from Tomo.ttime2D import TTime2D
from Station.station import Station
from Tomo.inverse import Inverse

from scipy.spatial import ConvexHull
from matplotlib import path

class Tomo:

    def __init__(self, params):

        self.params = params

        self.LAT0 = float(params['box'][0])
        self.LON0 = float(params['box'][1])
        self.RADIUS = float(params['box'][2])

        self.VM = Velmod(params['box'])

        self.st = Station(params['stafile'])

        # Dispersion curve variables
        self.STA1 = []
        self.STA2 = []
        self.P1 = []
        self.P2 = []
        self.PER = []
        self.TIME = []
        self.CMP =[]

        print("TOMO")
        print(params)
        # Read dispersion curves
        for root, dirnames, filenames in sorted(os.walk(self.params['dispdir'])):
            for filename in filenames:
                file = os.path.join(root, filename)
                #print("file=",file)
                if fnmatch.fnmatch(filename, 'disp_??_*_*_*.dat'):

                    fspl = filename.split('_')
                    self.CMP.append(fspl[1])
                    self.STA1.append(fspl[2])
                    self.STA2.append(fspl[3])
                    #print("Reading ",file)
                    f=open(file,'r')
                    print(file)
                    dist=float(f.readline())
                    f.close()
                    d=np.loadtxt(file, skiprows=1)
                    per=d[:,0]
                    self.PER.append(per)
                    self.P1.append(np.min(per))
                    self.P2.append(np.max(per))
                    vel=d[:,1]
                    self.TIME.append(dist/vel)
                    #print("Read ", filename, " N=",vel.size)

        ndat = len(self.STA1)

        self.X1 = np.empty(ndat)
        self.Y1 = np.empty(ndat)
        self.X2 = np.empty(ndat)
        self.Y2 = np.empty(ndat)

        # Stations in the box reference system
        for i in range(ndat):

            la, lo = self.st.getll(self.STA1[i])
            self.X1[i], self.Y1[i], zz, ll = utm.from_latlon(la, lo, self.VM.Z0)
            self.X1[i] = self.X1[i] - self.VM.X0
            self.Y1[i] = self.Y1[i] - self.VM.Y0
            #print(self.STA1[i], la, lo, self.X1[i], self.Y1[i])

            la, lo = self.st.getll(self.STA2[i])
            self.X2[i], self.Y2[i], zz, ll = utm.from_latlon(la, lo, self.VM.Z0)
            self.X2[i] = self.X2[i] - self.VM.X0
            self.Y2[i] = self.Y2[i] - self.VM.Y0
            #print(self.STA2[i], la, lo, self.X2[i], self.Y2[i])

    def get_times(self, comp, per):

        # Return observed travel times for a given period
        print("get_times ",comp,per)

        # Returns XS1, YS1, XS2, YS2, TS, avgv
        XS1 = []
        YS1 = []
        XS2 = []
        YS2 = []
        TS = []

        sd=0.
        st=0.
        for i in range(len(self.STA1)):
            if self.CMP[i]==comp and per>=self.P1[i] and per<=self.P2[i]:
                XS1.append(self.X1[i])
                XS2.append(self.X2[i])
                YS1.append(self.Y1[i])
                YS2.append(self.Y2[i])
                tt = np.interp(per,self.PER[i],self.TIME[i])
                TS.append(tt)
                dist = np.sqrt( (XS1[-1]-XS2[-1])**2 + (YS1[-1]-YS2[-1])**2 )
                sd = sd + dist
                st = st + tt

        if len(TS)<10:
            return XS1, YS1, XS2, YS2, TS, -1

        avgv = sd / st

        XS1 = np.array(XS1)
        XS2 = np.array(XS2)
        YS1 = np.array(YS1)
        YS2 = np.array(YS2)
        TS = np.array(TS)

        # Check for stations inside box
        size = self.VM.size
        id1 = (XS1>-size)
        id2 = np.multiply(id1, XS1<size)
        id3 = np.multiply(id2, YS1>-size)
        id4 = np.multiply(id3, YS1<size)

        id5 = np.multiply(id4, XS2>-size)
        id6 = np.multiply(id5, XS2<size)
        id7 = np.multiply(id6, YS2>-size)
        id8 = np.multiply(id7, YS2<size)

        idx = np.argwhere(id8)
        XS1 = np.squeeze(XS1[idx])
        XS2 = np.squeeze(XS2[idx])
        YS1 = np.squeeze(YS1[idx])
        YS2 = np.squeeze(YS2[idx])
        TS = np.squeeze(TS[idx])

        return XS1, YS1, XS2, YS2, TS, avgv

    def getrays(self, XS1, YS1, XS2, YS2, sc, m):

        NT = len(XS1)
        rays = [None] * NT

        if self.params['ttcalc']=='STRAIGHT':

            NINT=100
            for k in range(NT):
                XP = np.linspace(XS1[k], XS2[k], NINT)
                YP = np.linspace(YS1[k], YS2[k], NINT)
                rr = np.vstack((XP,YP)).T
                rays[k] = rr

        else:

            print("TTIME2D calcrays", flush=True)
            self.VM.coef2mod(sc, m)
            DL = (2*self.RADIUS)/self.params['ttgrid']
            TT = TTime2D(self.params['ttgrid'],self.params['ttgrid'], DL, -self.RADIUS, -self.RADIUS)
            TT.setvelmod(self.VM)

            count=0
            #import psutil
            #print(psutil.virtual_memory())
            flag=np.full(NT, dtype=bool, fill_value=True)
            while sum(flag)>0:

                # Search first empty source
                for i in range(NT):
                    if flag[i]: break

                xs = XS1[i]
                ys = YS1[i]
                TT.calc(xs,ys)
                for i in range(NT):
                    if XS1[i]==xs and YS1[i]==ys:

                        rx, ry, rt = TT.backray(XS2[i], YS2[i])
                        flag[i]=False
                        rr = np.vstack((rx,ry)).T
                        rays[i] = rr

                        count=count+1

        return rays

    def getsynth(self, NT, sc, m, rays):

        TS = np.empty(NT)
        S = self.VM.coef2mod(sc, m)

        for i in range(NT):
            rx = rays[i][:, 0]
            ry = rays[i][:, 1]
            rxr = np.roll(rx, -1)
            ryr = np.roll(ry, -1)
            dist = np.sqrt( (rx-rxr)**2 + (ry-ryr)**2 )

            slo = S(rx, ry)

            slor = np.roll(slo, -1)
            slom = (slo+slor)/2
            dt = dist*slom
            TS[i]=np.sum(dt[:-1])

        return TS

    def writeres(self, XS1, YS1, XS2, YS2, TS, comp, per, sc, mod, RSS0):

        rays = self.getrays(XS1, YS1, XS2, YS2, sc, mod)
        sy = self.getsynth(len(XS1), sc, mod, rays)

        # Compute residuals
        res = TS - sy
        rss = np.sum(res ** 2)

        # Write data, synth & residuals
        f=open("%s/res_%s_%3.2f_%d.dat" % (self.params['tomodir'],comp,per,sc), 'w')

        f.write("RSS %f %f %f\n" % (rss, RSS0, 100*(rss-RSS0)/RSS0))

        for i in range(len(XS1)):

            la1, lo1 = utm.to_latlon(XS1[i] + self.VM.X0, YS1[i] + self.VM.Y0, self.VM.Z0, self.VM.L0)
            la2, lo2 = utm.to_latlon(XS2[i] + self.VM.X0, YS2[i] + self.VM.Y0, self.VM.Z0, self.VM.L0)
            dist = np.sqrt( (XS1[i]-XS2[i])**2 + (YS1[i]-YS2[i])**2 )
            vobs = dist/TS[i]
            vcalc = dist/sy[i]
            f.write("%f %f %f %f  %f %f %f %f  %f  %f %f  %f %f\n" % (XS1[i],YS1[i],XS2[i],YS2[i], la1,lo1, la2,lo2, dist/1000., TS[i], sy[i], vobs/1000., vcalc/1000.))

        f.close()

    def writerays(self, XS1, YS1, XS2, YS2, comp, per, sc, mod):

        rays = self.getrays(XS1, YS1, XS2, YS2, sc, mod)

        # Write rays
        f = open("%s/rays_%s_%3.2f_%d.dat" % (self.params['tomodir'], comp, per, sc), 'w')
        for i in range(len(XS1)):
            rx = rays[i][:, 0]
            ry = rays[i][:, 1]
            for j in range(len(rx)):
                #print(rx[j],ry[j],rx[j]+self.VM.X0, ry[j]+self.VM.Y0)
                la, lo = utm.to_latlon(rx[j]+self.VM.X0, ry[j]+self.VM.Y0, self.VM.Z0, self.VM.L0)
                f.write("%f %f  %f %f\n" % (lo, la, rx[j], ry[j]))
            f.write(">\n")
        f.close()

    def writemod(self, ofile, DX, sc, mod, XS1, YS1, XS2, YS2):

        self.VM.coef2mod(sc, mod)

        xg, yg = self.VM.getnodes(sc)

        # Convex hull
        XP = np.hstack((XS1, XS2))
        YP = np.hstack((YS1, YS2))
        PT=np.vstack((XP,YP)).T

        hull = ConvexHull(PT)
        ch=hull.points[hull.vertices,:]

        xmin = np.min(ch[:, 0])
        xmax = np.max(ch[:, 0])
        ymin = np.min(ch[:, 1])
        ymax = np.max(ch[:, 1])

        pol=[]
        for i in range(ch.shape[0]):
            pol.append((ch[i,0],ch[i,1]))
        pol.append((ch[0, 0], ch[0, 1]))

        pa = path.Path(pol, closed=True)

        f=open(ofile, 'w')
        for xo in np.arange(-self.VM.size, self.VM.size, DX):
            for yo in np.arange(-self.VM.size, self.VM.size, DX):
                if xo<xmin or xo>xmax or yo<ymin or yo>ymax: continue
                if pa.contains_point((xo,yo)):
                    la, lo = utm.to_latlon(xo+self.VM.X0, yo+self.VM.Y0, self.VM.Z0, self.VM.L0)
                    vel = self.VM.getV(xo,yo)
                    f.write("%f %f  %f %f  %f\n" % (la, lo, xo, yo, vel/1000.) )
        f.close()

    def writenodes(self, ofile, sc):

        xg, yg = self.VM.getnodes(sc)

        f = open(ofile, 'w')
        for i in range(len(xg)):
            la, lo = utm.to_latlon(xg[i] + self.VM.X0, yg[i] + self.VM.Y0, self.VM.Z0, self.VM.L0)
            f.write("%f %f  %f %f\n" % (la, lo, xg[i], yg[i]))
        f.close()

    def scaleupgrade(self, m, sc1, sc2):

        print(">>>>>>> scaleupgrade")
        V = self.VM.coef2mod(sc1, m)
        m2 = self.VM.mod2coef(V, sc2)

        return m2

    def dotomo(self):

        Inv = Inverse(self.params)

        # Create tomodir
        if not os.path.isdir(self.params['tomodir']):
            os.mkdir(self.params['tomodir'])

        for cmp in self.params['comp']:
            for per in self.params['periods']:

                print("*********************************************")
                print("comp=",cmp," period=",per)

                # Observed data
                XS1, YS1, XS2, YS2, TS, avgv = self.get_times(cmp,per)

                N=len(TS)
                print(">>>>> N=",N," avgv=",avgv)
                #print(XS1,YS1,XS2,YS2)
                if N<10: continue

                # LINEAR
                if self.params['inversion']=="LINEAR":

                    print("LINEAR INVERSION - STRAIGHT-LINE RAY TRACING - NO MULTISCALE", flush=True)

                    sc=self.params['scale']
                    self.VM.setType('lin')
                    Inv.settickhonov(self.VM.gettikhonov(sc))

                    # Reference model
                    print("avgv=",avgv)
                    m0 = self.VM.homo2coef(sc, avgv)

                    rays = self.getrays(XS1, YS1, XS2, YS2, sc, m0)
                    T0 = self.getsynth(N, sc, m0, rays)
                    res = T0 - TS
                    RSS0 = np.sum(res ** 2)

                    # Kernel
                    DELTAM = 1e-6
                    M = m0.size
                    print("KERNEL N=", N, " M=", M, flush=True)
                    G = np.empty((N, M))
                    for k in range(M):
                        # Calc residuals and G
                        mod = deepcopy(m0)
                        mod[k] = mod[k] + DELTAM
                        tk = self.getsynth(N, sc, mod, rays)
                        tk = (tk - T0) / DELTAM
                        G[:, k] = tk.T

                    #print("G=",G.shape)
                    val = lambda m : self.VM.check_valid(sc,m)
                    dm, flag = Inv.getinverse(G, m0, res, val, 'LINEAR')
                    if not flag:
                        print("ERROR during INVERSION!")
                        continue
                    mfinal = m0 + dm

                    print(">>>>> WRITEMOD")
                    ofile = "%s/tomo_%s_%3.2f_%d.dat" % (self.params['tomodir'], cmp, per, sc)
                    self.writemod(ofile, self.RADIUS/100., sc, mfinal, XS1, YS1, XS2, YS2)

                    print(">>>>> WRITENODES")
                    ofile = "%s/nodes_%s_%3.2f_%d.dat" % (self.params['tomodir'], cmp, per, sc)
                    self.writenodes(ofile,sc)

                    print(">>>>> WRITERES")
                    self.writeres(XS1, YS1, XS2, YS2, TS, cmp, per, sc, mfinal, RSS0)

                    print(">>>>> WRITERAYS", flush=True)
                    self.writerays(XS1, YS1, XS2, YS2, cmp, per, sc, mfinal)

                # NON-LINEAR
                elif not 'multiscale' in self.params.keys():

                    # No multiscale

                    print("NON-LINEAR INVERSION - SHORTEST-PATH RAY TRACING - NO MULTISCALE", flush=True)

                    sc = self.params['scale']
                    Inv.settickhonov(self.VM.gettikhonov(sc))

                    # Reference model
                    print("avgv=", avgv)
                    m0 = self.VM.homo2coef(sc, avgv)

                    ###################################################33
                    # Preliminary linear inversion
                    print("Preliminary linear inversion")

                    self.VM.setType('lin')

                    self.params['ttcalc'] = 'STRAIGHT'
                    rays = self.getrays(XS1, YS1, XS2, YS2, sc, m0)
                    self.params['ttcalc'] = 'SP'
                    T0 = self.getsynth(N, sc, m0, rays)
                    res = T0 - TS
                    RSS0 = np.sum(res ** 2)

                    # Kernel
                    DELTAM = 1e-6
                    M = m0.size
                    print("KERNEL N=", N, " M=", M, flush=True)
                    G = np.empty((N, M))
                    for k in range(M):
                        # Calc residuals and G
                        mod = deepcopy(m0)
                        mod[k] = mod[k] + DELTAM
                        tk = self.getsynth(N, sc, mod, rays)
                        tk = (tk - T0) / DELTAM
                        G[:, k] = tk.T

                    #print("G=", G.shape)
                    val = lambda m: self.VM.check_valid(sc, m)
                    dm, flag = Inv.getinverse(G, m0, res, val, 'LINEAR')
                    if not flag:
                        print("ERROR during INVERSION!")
                        continue
                    mlin = m0 + dm

                    m0 = mlin
                    dm0 = dm

                    #print("Preliminary linear model=",m0)

                    #################################
                    # MAIN LOOP FOR NON-LINEAR
                    print("MAIN LOOP FOR NON-LINEAR")

                    # Pass from linear to non-linear
                    V = self.VM.coef2mod(sc, m0)
                    #print("LINEAR m0=",m0)
                    self.VM.setType('log')
                    m0 = self.VM.mod2coef(V, sc)
                    #print("NONLINEAR m0=", m0)

                    mref=deepcopy(m0)

                    #m0 = self.VM.homo2coef(sc, avgv)
                    #print("NONLINEAR m0=", m0)
                    #exit(1)

                    deltam=1e10
                    iter = 0
                    lam = 1.0

                    while deltam>self.params['invdm'] and iter<self.params['invmaxiter']:

                        print("ITER=",iter, flush=True)

                        #print("CALCRAYS", flush=True)
                        rays = self.getrays(XS1, YS1, XS2, YS2, sc, m0)
                        T0 = self.getsynth(N, sc, m0, rays)
                        res = T0 - TS
                        RSS0 = np.sum(res ** 2)

                        # Kernel
                        DELTAM = 1e-2
                        M = m0.size
                        print("KERNEL N=", N, " M=", M, flush=True)
                        G = np.empty((N, M))
                        for k in range(M):
                            # Calc residuals and G
                            #print("KERNEL ",k,M)
                            mod = deepcopy(m0)
                            mod[k] = mod[k] + DELTAM
                            tk = self.getsynth(N, sc, mod, rays)
                            tk = (tk - T0) / DELTAM
                            G[:, k] = tk.T

                        val = lambda m: self.VM.check_valid(sc, m)
                        dm, flag = Inv.getinverse(G, m0, res, val, 'NONLINEAR', mref=mref)
                        dm = dm * lam
                        lam = lam / 2.
                        if not flag:
                            print("ERROR during INVERSION!")
                            break
                        m1 = m0 + dm

                        #print("mref ", mref)
                        #print("m0 ", m0)
                        #print("dm ", dm)
                        #print("m1 ", m1)

                        #deltam=100*np.abs(dm) / np.abs(m0)
                        #deltam[np.isnan(deltam)] = 0.
                        #deltam[np.isinf(deltam)] = 0.
                        #print("deltam ", deltam)
                        #deltam = np.max(deltam)
                        #print("iter=",iter," deltam=",deltam,"% dm=",np.linalg.norm(dm)/len(dm), flush=True)
                        deltam = 100*np.linalg.norm(dm)/np.linalg.norm(m0)
                        print("iter=", iter, " deltam=", deltam, flush=True)

                        """
                        print(">>>>> WRITEMOD PART")
                        ofile = "PART/tomo_%s_%3.2f_%d_%02d.dat" % (cmp, per, sc,iter)
                        self.writemod(ofile, self.RADIUS / 100., sc, m1, XS1, YS1, XS2, YS2)
                        """

                        #print("dm/dm0 = ",100.*np.max(ddm),'%')
                        #dm0 = dm

                        iter = iter+1

                        m0 = m1

                    mfinal = m0

                    print(">>>>> WRITEMOD")
                    ofile = "%s/tomo_%s_%3.2f_%d.dat" % (self.params['tomodir'], cmp, per, sc)
                    self.writemod(ofile, self.RADIUS / 100., sc, mfinal, XS1, YS1, XS2, YS2)

                    print(">>>>> WRITENODES")
                    ofile = "%s/nodes_%s_%3.2f_%d.dat" % (self.params['tomodir'], cmp, per, sc)
                    self.writenodes(ofile, sc)

                    print(">>>>> WRITERES")
                    self.writeres(XS1, YS1, XS2, YS2, TS, cmp, per, sc, mfinal, RSS0)

                    print(">>>>> WRITERAYS", flush=True)
                    self.writerays(XS1, YS1, XS2, YS2, cmp, per, sc, mfinal)


                else:

                    # Multiscale
                    print("NON-LINEAR INVERSION - SHORTEST-PATH RAY TRACING - MULTISCALE", flush=True)

                    # Preliminary model for scale 0
                    self.VM.setType('log')
                    print("avgv=", avgv)
                    mprel = self.VM.homo2coef(0, avgv)

                    m0 = self.scaleupgrade(mprel, 0, 1)

                    print("mprel=",mprel," m0=",m0)

                    ###################################################33

                    for sc in range(1,self.params['multiscale']+1):

                        print("MULTISCALE INVERSION scale=",sc, flush=True)

                        Inv.settickhonov(self.VM.gettikhonov(sc))

                        rays = self.getrays(XS1, YS1, XS2, YS2, sc, m0)
                        T0 = self.getsynth(N, sc, m0, rays)
                        res = T0 - TS
                        RSS0 = np.sum(res ** 2)

                        #################################
                        # MAIN NON-LINEAR LOOP
                        print("MAIN LOOP FOR NON-LINEAR", flush=True)

                        deltam = 1e10
                        iter = 0
                        lam = 1
                        mref = deepcopy(m0)

                        while deltam > self.params['invdm'] and iter < self.params['invmaxiter']:

                            print("ITER=", iter)

                            rays = self.getrays(XS1, YS1, XS2, YS2, sc, m0)
                            T0 = self.getsynth(N, sc, m0, rays)
                            res = T0 - TS
                            RSS0 = np.sum(res ** 2)

                            # Kernel
                            DELTAM = 1e-2
                            M = m0.size
                            print("KERNEL N=", N, " M=", M, flush=True)
                            G = np.empty((N, M))
                            for k in range(M):
                                mod = deepcopy(m0)
                                mod[k] = mod[k] + DELTAM
                                tk = self.getsynth(N, sc, mod, rays)
                                tk = (tk - T0) / DELTAM
                                G[:, k] = tk.T

                            # Inversion
                            val = lambda m: self.VM.check_valid(sc, m)
                            dm, flag = Inv.getinverse(G, m0, res, val, 'NONLINEAR', mref=mref)
                            dm = dm * lam
                            lam = lam / 2.
                            if not flag:
                                print("ERROR during INVERSION!")
                                break
                            m1 = m0 + dm

                            deltam = 100 * np.linalg.norm(dm) / np.linalg.norm(m0)
                            print("iter=", iter, " deltam=", deltam, flush=True)

                            iter = iter + 1

                            m0 = m1

                        mfinal = m0

                        print(">>>>> WRITEMOD")
                        ofile = "%s/tomo_%s_%3.2f_%d.dat" % (self.params['tomodir'], cmp, per, sc)
                        self.writemod(ofile, self.RADIUS / 100., sc, mfinal, XS1, YS1, XS2, YS2)

                        print(">>>>> WRITENODES")
                        ofile = "%s/nodes_%s_%3.2f_%d.dat" % (self.params['tomodir'], cmp, per, sc)
                        self.writenodes(ofile, sc)

                        print(">>>>> WRITERES")
                        self.writeres(XS1, YS1, XS2, YS2, TS, cmp, per, sc, mfinal, RSS0)

                        print(">>>>> WRITERAYS", flush=True)
                        self.writerays(XS1, YS1, XS2, YS2, cmp, per, sc, mfinal)

                        if sc < self.params['multiscale']:
                            print("scale upgrade")
                            m0 = self.scaleupgrade(mfinal, sc, sc + 1)
