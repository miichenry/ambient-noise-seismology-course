
import obspy
import fnmatch
import os
import datetime
import math
import numpy
from copy import deepcopy

from Station.station import Station

class Finder:

    # General info
    stations=[]     # This is taken from the seismograms
    delta=[]
    nptw=-1

    # Station infos
    st=[]

    # File info
    folder=[]
    fname=[]
    fpos=[]
    fsta=[]
    fcmp=[]
    fstarttime=[]
    fendtime=[]
    fnpts=[]

    # Window info
    twin=[]     # Starting time of each window
    swin=[]     # Matrix of booleans, True=station has the windows
    # Matrices with index reference to files for each window
    ifilez = []
    ifilee = []
    ifilen = []

    def __init__(self, datadir, st):

        self.folder=datadir

        self.st=st

        for root, dirnames, filenames in os.walk(datadir):
            for filename in filenames:
                print(root,filename)
                file = os.path.join(root, filename)
                try:
                    st = obspy.read(file)
                    for i in range(len(st)):
                        tr = st[i]
                        if not tr.stats.station in self.st.staname: continue
                        self.fname.append(file)
                        self.fpos.append(i)
                        self.fsta.append(tr.stats.station)
                        self.fcmp.append(tr.stats.channel[2])
                        self.fstarttime.append(tr.stats.starttime)
                        self.fendtime.append(tr.stats.endtime)
                        self.fnpts.append(tr.stats.npts)
                        self.delta.append(tr.stats.delta)

                except:
                    pass

        self.stations = sorted(list(set(self.fsta)))
        self.delta = sorted(list(set(self.delta)))
        if len(self.delta)>1:
            print("ERROR: different sampling rates in input files!")
            raise Exception
        self.delta=self.delta[0]

    def check_windows(self, wlen_=200, wpwr_=-1):

        if wpwr_==-1:
            npts = numpy.floor(wlen_/self.delta)
            wpwr = numpy.ceil(numpy.log2(npts))
            self.nptw = 2 ** wpwr
            wlen = self.nptw * self.delta
            self.wlen = wlen

        else:
            self.nptw=2**wpwr_
            wlen=self.nptw*self.delta
            self.wlen=wlen

        t1=min(self.fstarttime)
        i1=self.fstarttime.index(t1)
        t2=max(self.fendtime)
        dt=t2-t1
        self.nwin=int(dt/wlen)

        self.twin=[0] * self.nwin
        for i in range(self.nwin):
            self.twin[i]=self.fstarttime[i1]+datetime.timedelta(seconds=i*wlen)

        swin_z = numpy.zeros((self.nwin, len(self.stations)), dtype=bool)
        swin_e = numpy.zeros((self.nwin, len(self.stations)), dtype=bool)
        swin_n = numpy.zeros((self.nwin, len(self.stations)), dtype=bool)

        self.ifilez = numpy.full((self.nwin, len(self.stations)), -1)
        self.ifilee = numpy.full((self.nwin, len(self.stations)), -1)
        self.ifilen = numpy.full((self.nwin, len(self.stations)), -1)

        for i in range(len(self.fsta)):

            ist=self.stations.index(self.fsta[i])

            i1 = int(math.ceil((self.fstarttime[i] - t1) / wlen))
            i2 = int(math.floor((self.fendtime[i] - t1) / wlen))

            if self.fcmp[i]=='Z':
                for j in range(i1, i2):
                    swin_z[j, ist] = True
                    self.ifilez[j, ist] = i
            elif self.fcmp[i]=='E':
                for j in range(i1, i2):
                    swin_e[j, ist] = True
                    self.ifilee[j, ist] = i
            elif self.fcmp[i]=='N':
                for j in range(i1, i2):
                    swin_n[j, ist] = True
                    self.ifilen[j, ist] = i
            else:
                print("ERROR: bad component "+self.fcmp[i]+" in file "+self.fname[i]+"!")
                raise Exception

        self.swin = numpy.multiply(swin_z, swin_e)
        self.swin = numpy.multiply(self.swin, swin_n)

    def _getwinz(self, win, ista):

        idx = self.ifilez[win, ista]
        st = obspy.read(self.fname[idx])
        tr = st[self.fpos[idx]]

        dt=self.twin[win]-self.fstarttime[idx]
        isam = int(dt / self.delta)

        return tr.data[isam:isam+self.nptw]

    def _getwine(self, win, ista):

        idx = self.ifilee[win, ista]
        st = obspy.read(self.fname[idx])
        tr = st[self.fpos[idx]]

        dt=self.twin[win]-self.fstarttime[idx]
        #isam=int(dt/(self.nptw*self.delta))
        isam = int(dt / self.delta)

        return tr.data[isam:isam+self.nptw]

    def _getwinn(self, win, ista):

        idx = self.ifilen[win, ista]
        st = obspy.read(self.fname[idx])
        tr = st[self.fpos[idx]]

        dt=self.twin[win]-self.fstarttime[idx]
        #isam=int(dt/(self.nptw*self.delta))
        isam = int(dt / self.delta)

        return tr.data[isam:isam+self.nptw]

    def getStreamZ(self, win):

        nsta=sum(self.swin[win,:])
        st=obspy.Stream()

        s=0
        for i in range(self.swin.shape[1]):
            if self.swin[win,i]:
                data=self._getwinz(win,i)
                s=s+1
                hdr=obspy.core.trace.Stats()
                hdr.delta=self.delta
                hdr.npts=len(data)
                hdr.starttime=self.twin[win]
                #hdr.endtime=self.twin[win]+hdr.delta*hdr.npts
                st=st+obspy.Trace(data, hdr)

        return st

    def getnwindows(self, sta1, sta2):

        i1 = self.stations.index(sta1)
        i2 = self.stations.index(sta2)

        count=0
        for i in range(len(self.twin)):
            if self.swin[i,i1] and self.swin[i,i2]:
                count=count+1

        return count

    def gettwindow(self, sta1, sta2, k):

        i1 = self.stations.index(sta1)
        i2 = self.stations.index(sta2)

        count=0
        for i in range(len(self.twin)):
            if self.swin[i, i1] and self.swin[i, i2]:
                if count==k:
                    break
                else:
                    count = count + 1

        return self.twin[i]

    def getwindow(self, sta1, sta2, k):

        i1 = self.stations.index(sta1)
        i2 = self.stations.index(sta2)

        count = 0
        for i in range(len(self.twin)):
            if self.swin[i, i1] and self.swin[i, i2]:
                if count==k:
                    break
                else:
                    count = count + 1

        w1z = self._getwinz(i, i1)
        w2z = self._getwinz(i, i2)
        w1e = self._getwine(i, i1)
        w2e = self._getwine(i, i2)
        w1n = self._getwinn(i, i1)
        w2n = self._getwinn(i, i2)

        # Rotate horizontal components
        lat1 = self.st.stalat[i1]
        lon1 = self.st.stalon[i1]
        lat2 = self.st.stalat[i2]
        lon2 = self.st.stalon[i2]

        dl=lon2-lon1
        X=math.cos(lat2*math.pi/180.)*math.sin(dl*math.pi/180.)
        Y=math.cos(lat1*math.pi/180.)*math.sin(lat2*math.pi/180.) - math.sin(lat1*math.pi/180.)*math.cos(lat2*math.pi/180.)*math.cos(dl*math.pi/180.)
        az=math.atan2(X, Y)

        w1r = w1n * math.cos(az) + w1e * math.sin(az)
        w2r = w2n * math.cos(az) + w2e * math.sin(az)

        w1t = w1n * math.sin(az) - w1e * math.cos(az)
        w2t = w2n * math.sin(az) - w2e * math.cos(az)

        return [w1z, w1r, w1t, w2z, w2r, w2t]

    def reportwin(self, outfile):

        wlen=self.nptw*self.delta

        f=open(outfile,'w')
        for i in range(len(self.twin)):
            f.write('%s %s\n' % (self.twin[i], self.twin[i]+datetime.timedelta(seconds=wlen)) )
            for j in range(len(self.stations)):
                if self.swin[i,j]:
                    iz = self.ifilez[i,j]
                    f.write('%s Z %s %s\n' % (self.fsta[iz], self.fname[iz], self.fpos[iz]) )
                    ie = self.ifilee[i, j]
                    f.write('%s E %s %s\n' % (self.fsta[ie], self.fname[ie], self.fpos[ie]))
                    i_n = self.ifilez[i, j]
                    f.write('%s N %s %s\n' % (self.fsta[i_n], self.fname[i_n], self.fpos[i_n]))

        f.close()
