
import numpy as np
from scipy.interpolate import griddata
from copy import deepcopy

from TTime.Velmod import Velmod2D

import matplotlib.pyplot as plt

class TTime2D:

    def __init__(self, NX, NZ, DL, X0=0, Z0=0):

        # Grid
        self.NX=NX
        self.NZ=NZ
        self.DL=DL

        self.X0=X0
        self.Z0=Z0

        # Point distribution within a cell
        self.XP = np.array([0,      0,          DL,     DL,         DL / 3, 2 * DL / 3, DL / 3, 2 * DL / 3])
        self.ZP = np.array([DL / 3, 2 * DL / 3, DL / 3, 2 * DL / 3, 0,      0,          DL,     DL])

        # Distance matrix
        DM = np.empty((8, 8))
        for i in range(8):
            DM[i, :] = np.sqrt((self.XP - self.XP[i]) ** 2 + (self.ZP - self.ZP[i]) ** 2)

        # Point distribution in auxiliary cell
        self.AXP = np.array([0,     -DL/3, 0,    DL/3])
        self.AZP = np.array([-DL/3, 0,     DL/3, 0])

        # Traveltime on grid points
        self.NPT = (NX + 1) * (2 * NZ) + (2 * NX) * (NZ + 1)
        self.PX = np.empty(self.NPT) # X coordinate
        self.PZ = np.empty(self.NPT) # Z coordinate

        #print("NPT=",self.NPT)

        ########################
        # Create grid structure
        ########################

        self.GS = np.empty((NX, NZ, 8), dtype=np.int32)

        # Loop over cells
        k = 0
        for i in range(NX):

            Xi = self.X0 + i*self.DL

            for j in range(NZ):

                #print("i,j,k=",i,j,k)

                Zi = self.Z0 + j * self.DL

                if i == 0:
                    # Left 0 (only first column)
                    self.PX[k] = Xi + self.XP[0]
                    self.PZ[k] = Zi + self.ZP[0]

                    self.GS[i, j, 0] = k
                    k = k + 1

                    # Left 1 (only first column)
                    self.PX[k] = Xi + self.XP[1]
                    self.PZ[k] = Zi + self.ZP[1]
                    self.GS[i, j, 1] = k
                    k = k + 1

                # Right 2
                self.PX[k] = Xi + self.XP[2]
                self.PZ[k] = Zi + self.ZP[2]
                self.GS[i, j, 2] = k
                if i < NX - 1: self.GS[i + 1, j, 0] = k
                k = k + 1

                # Right 3
                self.PX[k] = Xi + self.XP[3]
                self.PZ[k] = Zi + self.ZP[3]
                self.GS[i, j, 3] = k
                if i < NX - 1: self.GS[i + 1, j, 1] = k
                k = k + 1

                if j==0:

                    # Bottom 4 (only first row)
                    self.PX[k] = Xi + self.XP[4]
                    self.PZ[k] = Zi + self.ZP[4]
                    self.GS[i, j, 4] = k
                    k = k + 1

                    # Bottom 5 (only first row)
                    self.PX[k] = Xi + self.XP[5]
                    self.PZ[k] = Zi + self.ZP[5]
                    self.GS[i, j, 5] = k
                    k = k + 1

                # Top 6
                self.PX[k] = Xi + self.XP[6]
                self.PZ[k] = Zi + self.ZP[6]
                self.GS[i, j, 6] = k
                if j < NZ - 1: self.GS[i, j + 1, 4] = k
                k = k + 1

                # Top 7
                self.PX[k] = Xi + self.XP[7]
                self.PZ[k] = Zi + self.ZP[7]
                self.GS[i, j, 7] = k
                if j < NZ - 1: self.GS[i, j + 1, 5] = k
                k = k + 1

        """
        # Print struct
        print("MAIN CELLS")
        for i in range(NX):
            for j in range(NZ):
                print(">>>>> i,j=",i, j, " idx=",self.GS[i,j,:])
                for l in range(8):
                    k=self.GS[i,j,l]
                    print("l=",l," k=",k," PX=",self.PX[k]," PZ=",self.PZ[k])
        """

        ########################
        # Create forward star
        ########################

        # Forward star
        self.FS = []

        # Fill FS and FSD
        for i in range(self.NPT):
            self.FS.append([])

        # Loop over cells
        for i in range(NX):
            for j in range(NZ):

                # Link within the cell i,j

                for a in range(8):
                    k = self.GS[i,j,a]
                    for b in range(8):
                        # Link a with b except for adjacent nodes
                        if b!=a:# and b !=forbidden[a]:
                            #self.FS[k].append((self.GS[i, j, b],i,j,a,b))
                            # 0 idx, 1 i cell, 2 j cell, 3 dist
                            self.FS[k].append((self.GS[i, j, b],i,j,DM[a,b]))

        """
        # Print FS
        print("FORWARD STAR")
        for i in range(self.NPT):
            print(i,self.PX[i],self.PZ[i])
            for j in range(len(self.FS[i])):
                f=self.FS[i][j]
                print("   ",j,f,self.PX[f[0]],self.PZ[f[0]])
        """

        # Create auxiliary grid structure
        self.AS = np.empty((NX-1, NZ-1, 4), dtype=np.int32)
        k = 0
        for i in range(NX - 1):
            for j in range(NZ - 1):

                self.AS[i, j, 0] = self.GS[i, j, 3]
                self.AS[i, j, 1] = self.GS[i, j, 7]
                self.AS[i, j, 2] = self.GS[i+1, j+1, 0]
                self.AS[i, j, 3] = self.GS[i+1, j+1, 4]

        """
        print("AUX CELLS")
        for i in range(NX-1):
            for j in range(NZ-1):
                print(">>>>> i,j=", i, j, " idx=", self.AS[i, j, :])
                for l in range(4):
                    k = self.AS[i, j, l]
                    print("l=", l, " k=", k, " PX=", self.PX[k], " PZ=", self.PZ[k])
        """

    def setvelmod(self, vm, type):

        self.type = type

        # Velocity grid
        self.vel = np.empty((self.NX, self.NZ))

        for i in range(self.NX):
            xp = self.X0 + (i + 0.5) * self.DL
            for j in range(self.NZ):
                zp = self.Z0 + (j + 0.5) * self.DL
                if type=='P':
                    self.vel[i, j] = vm.getVP(xp, zp)
                else:
                    self.vel[i, j] = vm.getVS(xp, zp)

        # Weight of the forward star
        self.dtmin=1e80
        self.dtmax=0.
        self.FST = []
        for i in range(self.NPT):
            self.FST.append([])
            n = len(self.FS[i])
            for j in range(n):
                # In FS: 0 idx, 1 i cell, 2 j cell, 3 dist
                d = self.FS[i][j][3]
                v = self.vel[self.FS[i][j][1],self.FS[i][j][2]]
                dt = d/v
                self.FST[i].append(dt)
                if dt<self.dtmin: self.dtmin=dt
                if dt>self.dtmax: self.dtmax=dt

    def getcell(self, RX, RZ):

        fx = (RX - self.X0) / self.DL
        fz = (RZ - self.Z0) / self.DL

        ix = int(np.floor(fx))
        iz = int(np.floor(fz))

        # fx, fz = relative normalized position within the cell
        fx = (fx - ix)
        fz = (fz - iz)

        flagAux = False
        if fz < 1./3. - fx:
            # Lower left
            ai = ix - 1
            aj = iz - 1
            flagAux = True

        elif fz < -2./3. + fx:
            # Lower right
            ai = ix
            aj = iz - 1
            flagAux = True

        elif fz > 2./3. + fx:
            # Upper left
            ai = ix - 1
            aj = iz
            flagAux = True

        elif fz > 5./3. - fx:
            # Upper right
            ai = ix
            aj = iz
            flagAux = True

        else:
            ai = ix
            aj = iz

        if flagAux:
            fx = (RX - self.X0) / self.DL - (ai + 1)
            fz = (RZ - self.Z0) / self.DL - (aj + 1)

        return ai, aj, fx * self.DL, fz * self.DL, flagAux

    def calc(self, SX, SZ):

        self.SX = SX
        self.SZ = SZ

        self.TT = np.full(self.NPT, 1e80)                   # Traveltimes
        self.SP = np.full(self.NPT, -2, dtype=np.int32)     # Shortest-path graph

        # Node types: A fixed, B wavefront, C virgin
        NT = np.char.chararray(self.NPT)
        NT.fill(b'C')

        # Seed source
        self.isx, self.isz, self.sfx, self.sfz, self.sflagAux = self.getcell(SX, SZ)

        if not self.sflagAux:
            dist = np.sqrt((self.sfx - self.XP) ** 2 + (self.sfz - self.ZP) ** 2)
            k = np.squeeze(self.GS[self.isx, self.isz, :])
            for l in range(8):
                self.TT[k[l]] = dist[l] / self.vel[self.isx, self.isz]
                NT[k[l]] = b'B'
                self.SP[k[l]] = -1

        else:
            dist = np.sqrt((self.sfx - self.AXP) ** 2 + (self.sfz - self.AZP) ** 2)
            k = np.squeeze(self.AS[self.isx, self.isz, :])
            vel = (self.vel[self.isx, self.isz] + self.vel[self.isx+1, self.isz] + self.vel[self.isx, self.isz+1] + self.vel[self.isx+1, self.isz+1]) / 4
            for l in range(4):
                self.TT[k[l]] = dist[l] / vel
                NT[k[l]] = b'B'
                self.SP[k[l]] = -1

        # Nodes in each bucket
        bucksi = []
        nbucks = int(self.dtmax/self.dtmin) + 2
        for i in range(nbucks): bucksi.append([])

        # Bucket for each node
        bucks = np.full((self.NPT,), fill_value=None, dtype=object)

        # Init buckets
        tmin = 0
        idx = np.where(NT == b'B')[0]
        for i in range(len(idx)):
            tt = self.TT[idx[i]]
            it = int((tt-tmin)/self.dtmin)
            bucksi[it].append(idx[i])
            bucks[idx[i]] = bucksi[it]

        # Dijkstra algorithm
        iter = 0
        while True:

            #print("ITER ",iter,bucksi)

            # Check if there are any nodes left in buckets
            flag=False
            for i in range(nbucks):
                if len(bucksi[i])>0:
                    flag=True
                    break

            if not flag: break

            # Propagate from the nodes of the first bucket
            for k in range(len(bucksi[0])):

                ii = bucksi[0][k]
                nf = len(self.FS[ii])
                for i in range(nf):

                    fs = self.FS[ii][i][0]

                    #print("Prop from node ",ii," to node ",fs," with type ",NT[fs])

                    if NT[fs] == b'C':

                        # Update virgin node
                        NT[fs] = b'B'
                        dt = self.FST[ii][i]
                        nt = self.TT[ii] + dt
                        self.TT[fs] = nt
                        self.SP[fs] = ii

                        # Update buckets
                        it = int((nt-tmin)/self.dtmin)
                        #print("Prop from node ", ii, " to node ", fs, " with type ", NT[fs])
                        #print("tmin=",tmin," nt=",nt," it=",it," nbucks=",nbucks)
                        bucksi[it].append(fs)
                        bucks[fs] = bucksi[it]

                    elif NT[fs] == b'B':

                        dt = self.FST[ii][i]
                        nt = self.TT[ii] + dt
                        if nt<self.TT[fs]:

                            #print("Update B node")
                            it = int((self.TT[fs]-tmin)/self.dtmin)

                            # Update node
                            self.TT[fs] = nt
                            self.SP[fs] = ii

                            # Update buckets
                            nit = int((nt-tmin)/self.dtmin)
                            if bucksi[nit]!=bucks[fs]:
                                #print("it,nit ",it,nit)
                                #print(bucksi[it])
                                #print(bucksi[nit])
                                # Remove from the old bucket
                                bucksi[it].remove(fs)
                                # Add to the new one
                                bucksi[nit].append(fs)
                                # Update node
                                bucks[fs]=bucksi[nit]

                        #else: print("No prop")

                # 3) Move to A set
                NT[ii] = b'A'

            # Update buckets
            del bucksi[0]       # Remove first
            bucksi.append([])   # Add next one
            tmin=tmin+self.dtmin
            iter = iter + 1

    def gettime(self, RX, RZ):

        ai, aj, fx, fz, flagAux = self.getcell(RX, RZ)

        if ai==self.isx and aj==self.isz and flagAux==self.sflagAux:
            # Inside source cell
            if not flagAux:
                vel=self.vel[ai,aj]
            else:
                vel = (self.vel[ai, aj] + self.vel[ai+1, aj] + self.vel[ai, aj+1] + self.vel[ai+1, aj+1]) / 4

            dist=np.sqrt((RX-self.SX)**2 + (RZ-self.SZ)**2)
            return dist/vel

        if not flagAux:
            # MAIN CELL
            dist = np.sqrt((fx - self.XP) ** 2 + (fz - self.ZP) ** 2)
            idx = self.GS[ai, aj, :]
            tm = self.TT[idx]
            if np.min(dist)>0:
                num=np.dot(1./dist,tm)
                den=np.sum(1./dist)
                return num/den
            else:
                imin = np.argmin(dist)
                idx = self.GS[ai, aj, imin]
                return self.TT[idx]

        else:
            # AUXILIARY CELL
            dist = np.sqrt((fx - self.AXP) ** 2 + (fz - self.AZP) ** 2)
            idx = self.AS[ai, aj, :]
            tm = self.TT[idx]
            if np.min(dist) > 0:
                num = np.dot(1. / dist, tm)
                den = np.sum(1. / dist)
                return num / den
            else:
                imin = np.argmin(dist)
                idx = self.AS[ai, aj, imin]
                return self.TT[idx]

    def plot_struct(self, ray=False):

        import matplotlib.pyplot as plt

        plt.figure()

        if ray:
            for k in range(self.NPT):
                if self.SP[k]>=0:
                    plt.plot([self.PX[k],self.PX[self.SP[k]]],[self.PZ[k],self.PZ[self.SP[k]]],'k')
                elif self.SP[k]==-1:
                    plt.plot([self.PX[k], self.SX], [self.PZ[k], self.SZ], 'k')

        for k in range(self.NPT):
            if self.TT[k]<1e80:
                plt.plot(self.PX[k], self.PZ[k], 'k.')
            else:
                plt.plot(self.PX[k], self.PZ[k], 'kx')

        x1 = self.DL / 3
        x2 = self.NX*self.DL - x1
        y1 = self.DL / 3
        y2 = self.NX * self.DL - y1

        plt.plot(self.SX, self.SZ, 'bo')
        plt.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],'k:')
        plt.xlim([-self.DL, (self.NX+1)*self.DL])
        plt.ylim([-self.DL, (self.NZ+1)*self.DL])
        plt.show()

    def backray(self, RX0, RZ0):

        XR=[]
        ZR=[]
        TR=[]

        # Initial traveltime
        XR.append(RX0)
        ZR.append(RZ0)

        # Initial cell
        irx, irz, fx, fz, flagAux = self.getcell(RX0, RZ0)

        # Starts already in the source cell
        if irx==self.isx and irz==self.isz and flagAux==self.sflagAux:

            if not flagAux:
                vel = self.vel[irx,irz]
            else:
                vel = (self.vel[irx, irz] + self.vel[irx+1, irz] + self.vel[irx, irz+1] + self.vel[irx+1, irz+1]) / 4

            dist = np.sqrt( (RX0-self.SX)**2 + (RZ0-self.SZ)**2 )
            TR.append(dist/vel)
            XR.append(self.SX)
            ZR.append(self.SZ)
            TR.append(0.)
            return XR, ZR, TR

        # Initial cell
        T0=self.gettime(RX0, RZ0)
        TR.append(T0)

        # Compute local gradient in the cell
        if not flagAux:
            k = self.GS[irx, irz, :]
            G=np.empty((8,3))
            d=np.empty(8)
            for l in range(8):
                G[l, 0] = 1.
                G[l, 1] = self.XP[l] / self.DL
                G[l, 2] = self.ZP[l] / self.DL
                d[l] = self.TT[k[l]]

        else:
            k = self.AS[irx, irz, :]
            G=np.empty((4,3))
            d=np.empty(4)
            for l in range(4):
                G[l, 0] = 1.
                G[l, 1] = self.XP[l] / self.DL
                G[l, 2] = self.ZP[l] / self.DL
                d[l] = self.TT[k[l]]

        Gig = np.linalg.pinv(G)
        m = np.dot(Gig, d)
        grad = m[1:3]
        grad = -grad / np.linalg.norm(grad)

        # Search the point along the gradient
        lmax=0
        vmax=0
        for l in range(len(k)):
            v=-np.array([self.PX[k[l]]-RX0, self.PZ[k[l]]-RZ0])
            v=v/np.linalg.norm(v)
            vg=np.dot(v,grad)
            if vg<vmax:
                vmax=vg
                lmax=l

        kmin = k[lmax]

        XR.append(self.PX[kmin])
        ZR.append(self.PZ[kmin])
        TR.append(self.TT[kmin])

        pre = self.SP[kmin]
        while pre!=-1:

            kmin=pre
            XR.append(self.PX[kmin])
            ZR.append(self.PZ[kmin])
            TR.append(self.TT[kmin])
            pre = self.SP[kmin]

            if pre==-1:
                XR.append(self.SX)
                ZR.append(self.SZ)
                TR.append(0.)

        return XR, ZR, TR
