
import numpy as np
from scipy.interpolate import griddata
from copy import deepcopy

from Tomo.velmod import Velmod

#import matplotlib.pyplot as plt

class TTime2D:

    def __init__(self, NX, NZ, DL, X0=0, Z0=0, tol=1e-2):

        # Grid
        self.NX=NX
        self.NZ=NZ
        self.DL=DL

        self.X0=X0
        self.Z0=Z0

        self.tol = tol

        # Point distribution within a cell
        self.XP = np.array([0,      0,          DL,     DL,         DL / 3, 2 * DL / 3, DL / 3, 2 * DL / 3])
        self.ZP = np.array([DL / 3, 2 * DL / 3, DL / 3, 2 * DL / 3, 0,      0,          DL,     DL])

        # Distance matrix
        self.DM = np.empty((8, 8))
        for i in range(8):
            self.DM[i, :] = np.sqrt((self.XP - self.XP[i]) ** 2 + (self.ZP - self.ZP[i]) ** 2)

        # Point distribution in auxiliary cell
        self.AXP = np.array([0,     -DL/3, 0,    DL/3])
        self.AZP = np.array([-DL/3, 0,     DL/3, 0])

        # Traveltime on grid points
        self.NPT = (NX + 1) * (2 * NZ) + (2 * NX) * (NZ + 1)
        self.PX = np.empty(self.NPT) # X coordinate
        self.PZ = np.empty(self.NPT) # Z coordinate

        # Create grid structure
        self.GS = np.empty( (NX, NZ, 8), dtype=np.int32)
        k = 0
        for i in range(NX):

            Xi = self.X0 + i*self.DL

            for j in range(NZ):

                Zi = self.Z0 + j * self.DL

                # Left
                if i==0:
                    self.PX[k] = Xi + self.XP[0]
                    self.PZ[k] = Zi + self.ZP[0]
                    self.GS[i, j, 0] = k
                    k = k + 1
                    self.PX[k] = Xi + self.XP[1]
                    self.PZ[k] = Zi + self.ZP[1]
                    self.GS[i, j, 1] = k
                    k = k + 1

                # Right
                self.PX[k] = Xi + self.XP[2]
                self.PZ[k] = Zi + self.ZP[2]
                self.GS[i, j, 2] = k
                if i < NX - 1: self.GS[i + 1, j, 0] = k
                k = k + 1
                self.PX[k] = Xi + self.XP[3]
                self.PZ[k] = Zi + self.ZP[3]
                self.GS[i, j, 3] = k
                if i < NX - 1: self.GS[i + 1, j, 1] = k
                k = k + 1

                # Bottom
                if j==0:
                    self.PX[k] = Xi + self.XP[4]
                    self.PZ[k] = Zi + self.ZP[4]
                    self.GS[i, j, 4] = k
                    k = k + 1
                    self.PX[k] = Xi + self.XP[5]
                    self.PZ[k] = Zi + self.ZP[5]
                    self.GS[i, j, 5] = k
                    k = k + 1

                # Top
                self.PX[k] = Xi + self.XP[6]
                self.PZ[k] = Zi + self.ZP[6]
                self.GS[i, j, 6] = k
                if j < NZ - 1: self.GS[i, j + 1, 4] = k
                k = k + 1
                self.PX[k] = Xi + self.XP[7]
                self.PZ[k] = Zi + self.ZP[7]
                self.GS[i, j, 7] = k
                if j < NZ - 1: self.GS[i, j + 1, 5] = k
                k = k + 1

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
        # Print struct
        print("MAIN CELLS")
        for i in range(NX):
            for j in range(NZ):
                print(">>>>> i,j=",i, j, " idx=",self.GS[i,j,:])
                for l in range(8):
                    k=self.GS[i,j,l]
                    print("l=",l," k=",k," PX=",self.PX[k]," PZ=",self.PZ[k])

        print("AUX CELLS")
        for i in range(NX-1):
            for j in range(NZ-1):
                print(">>>>> i,j=", i, j, " idx=", self.AS[i, j, :])
                for l in range(4):
                    k = self.AS[i, j, l]
                    print("l=", l, " k=", k, " PX=", self.PX[k], " PZ=", self.PZ[k])
        """

    def setvelmod(self, vm):

        # Velocity grid
        self.vel = np.empty((self.NX, self.NZ))

        for i in range(self.NX):
            xp = self.X0 + (i + 0.5) * self.DL
            for j in range(self.NZ):
                zp = self.Z0 + (j + 0.5) * self.DL
                self.vel[i, j] = vm.getV(xp, zp)

        """
        import matplotlib.pyplot as plt
        plt.pcolor(self.vel.T)
        plt.colorbar()
        plt.show()
        """

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

    def updateCell(self, i, j, vrb=False):

        if i<0 or j<0 or i>=self.NX or j>=self.NZ:
            return False

        flag = False

        k = np.squeeze(self.GS[i, j, :])

        # Initial
        Ti = self.TT[k]

        # Compute traveltime matrix and shortest path
        TM = self.DM / self.vel[i, j]
        for l in range(8):
            TM[l, :] = TM[l, :] + Ti[l]

        Tf = Ti
        for l in range(8):
            imin = np.argmin(TM[:, l])
            NT = TM[imin, l]
            if NT < Tf[l]:
                if vrb: print("UPDATE ",i,j,l,np.abs(NT - Tf[l]))
                #if np.abs(2*(NT - Tf[l])/(NT + Tf[l]))>self.tol: flag=True
                if np.abs(NT - Tf[l]) > self.tol: flag = True
                #flag = True
                self.TT[k[l]] = NT
                self.SP[k[l]] = k[imin]

        return flag

    def calc(self, SX, SZ, maxiter=1):

        self.SX = SX
        self.SZ = SZ

        #print("CALC")
        #print(self.X0,self.Z0,self.NX,self.NZ,self.DL)
        #print(SX,SZ)
        #exit(1)

        self.TT = np.full(self.NPT, 1e10)                   # Traveltimes
        self.SP = np.full(self.NPT, -1, dtype=np.int32)     # Shortest-path graph

        # Seed source
        self.isx, self.isz, self.sfx, self.sfz, self.sflagAux = self.getcell(SX, SZ)

        if not self.sflagAux:
            dist = np.sqrt((self.sfx - self.XP) ** 2 + (self.sfz - self.ZP) ** 2)
            k = np.squeeze(self.GS[self.isx, self.isz, :])
            for l in range(8):
                self.TT[k[l]] = dist[l] / self.vel[self.isx, self.isz]

        else:
            dist = np.sqrt((self.sfx - self.AXP) ** 2 + (self.sfz - self.AZP) ** 2)
            k = np.squeeze(self.AS[self.isx, self.isz, :])
            vel = (self.vel[self.isx, self.isz] + self.vel[self.isx+1, self.isz] + self.vel[self.isx, self.isz+1] + self.vel[self.isx+1, self.isz+1]) / 4
            for l in range(4):
                self.TT[k[l]] = dist[l] / vel

        rmax = int( max(self.isx, self.isz, self.NX - self.isx - 1, self.NZ - self.isz - 1) * 1.5 )
        cx = (SX - self.X0) / self.DL
        cz = (SZ - self.Z0) / self.DL

        flag = True
        t = 0
        while flag and t<maxiter:

            flag = False

            # Expanding circle (clockwise)
            for r in range(0, rmax):
                nang = max(r * 64, 1)
                dang = 2 * np.pi / nang

                px = -1
                pz = -1

                for k in range(nang):
                    ix = int(cx + r * np.cos(k * dang))
                    iz = int(cz + r * np.sin(k * dang))

                    if ix != px or iz != pz:
                        f1 = self.updateCell(ix, iz)
                        flag = flag or f1

                    px = ix
                    pz = iz

            # Expanding circle (anticlockwise)
            for r in range(0, rmax):
                nang = max(r * 64, 1)
                dang = 2 * np.pi / nang

                px = -1
                pz = -1

                for k in range(nang, 0, -1):
                    ix = int(cx + r * np.cos(k * dang))
                    iz = int(cz + r * np.sin(k * dang))

                    if ix != px or iz != pz:
                        f1 = self.updateCell(ix, iz)
                        flag = flag or f1

                    px = ix
                    pz = iz

            t = t + 1

        #print("t=",t)

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

        if ray:
            for k in range(self.NPT):
                if self.SP[k]!=-1:
                    plt.plot([self.PX[k],self.PX[self.SP[k]]],[self.PZ[k],self.PZ[self.SP[k]]],'k')
                else:
                    plt.plot([self.PX[k], self.SX], [self.PZ[k], self.SZ], 'k')

        for k in range(self.NPT):
            if self.TT[k]<1e10:
                pass
                #plt.plot(self.PX[k], self.PZ[k], 'ko')
            else:
                plt.plot(self.PX[k], self.PZ[k], 'ro')

        plt.plot(self.SX, self.SZ, 'bo')
        plt.xlim([0, self.NX])
        plt.ylim([0, self.NZ])
        #plt.show()

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
