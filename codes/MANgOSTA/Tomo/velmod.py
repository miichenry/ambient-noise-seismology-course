
import numpy as np
import utm
from scipy.interpolate import LinearNDInterpolator

from scipy.spatial import ConvexHull
from matplotlib import path

class Velmod:

# Performs linear interpolation of the slowness between model nodes

# Every time the model is updated its interpolation function S is updated as well

# def getVel(self, x, y): returns velocioty (through linear interpolation of slowness)

# Type='lin' for coefficients are in slowness
# Type='log' for coefficients are in logarithm of slowness

# def checkinside(self, X, Y):

# def getncoefsc(self, sc):

# def homo2coef(self, sc, val):

# def mod2coef(self, S, sc):
# returns coefficients from an arbitrary interpolation function S of the slowness

# def coef2mod(self, sc, m):
# returns an interpolation function of the slowness from coefficients

    def __init__(self, box, type='lin'):

        self.LAT0 = float(box[0])
        self.LON0 = float(box[1])
        self.X0, self.Y0, self.Z0, self.L0 = utm.from_latlon(self.LAT0, self.LON0)
        self.size=float(box[2])

        # Define boundaries
        self.XMIN = self.X0 - self.size
        self.XMAX = self.X0 + self.size
        self.YMIN = self.Y0 - self.size
        self.YMAX = self.Y0 + self.size

        self.type=type
        #self.wav=wav

        #SCMAX=7
        #self.DW=DWTO2(SCMAX, wav)

    def setType(self, type):

        self.type=type

    def checkinside(self, X, Y):
        if X < self.XMIN or X > self.XMAX or Y < self.YMIN or Y > self.YMAX:
            return False
        return True

    def getncoefsc(self, sc):

        if self.type=='lin' or self.type=='log':
            if sc==0: N=1
            else:
                N = (2**(sc-1)+1)**2

        return N

        #elif self.type == 'wav':
        #    return self.DW.gettotcoefscale(sc)

    def homo2coef(self, sc, avgv):

        #print("homo2coef type=",self.type," avgv=",avgv)

        if self.type=='lin':
            if sc==0:
                self.coefs = np.array([1./avgv])
            else:
                N = 2**(sc-1) + 1
                grid=np.full((N,N), fill_value=1./avgv)
                self.coefs = np.reshape(grid,grid.size)

        elif self.type=='log':
            ll=np.log(1./avgv)
            if sc==0:
                self.coefs = np.array([ll])
            else:
                N = 2**(sc-1) + 1
                grid=np.full((N,N), fill_value=ll)
                self.coefs = np.reshape(grid,grid.size)

        #print("self.coefs=",self.coefs)

        """
        elif self.type=='wav':
            self.DW.sethomoval(val)
            #print("COEF=",self.DW.coef)
            ntcoef=self.DW.gettotcoefscale(sc)
            #print("sc=",sc," ntcoef=",ntcoef)
            m=np.empty(ntcoef)
            k=0
            for i in range(sc+1):
                ncoef=self.DW.getncoefscale(sc)
                #print("i=",i,ncoef,k)
                for j in range(ncoef):
                    m[k]=self.DW.getcoefsc(sc,j)
                    #print("j=",j,k)
                    k=k+1
            return m
        """

        self.coef2mod(sc, self.coefs)

        return self.coefs

    def mod2coef(self, S, sc):
        # S is an interpolator of the slowness

        if self.type=='lin':
            if sc == 0:
                self.coefs = np.array([S(0,0)])

            else:
                N = 2 ** (sc - 1) + 1
                ag = np.linspace(-self.size, self.size, N)
                xg, yg = np.meshgrid(ag, ag)
                xg = np.reshape(xg, xg.size)
                yg = np.reshape(yg, yg.size)
                self.coefs = S(xg, yg)

            #print("LIN self.coefs=",self.coefs)

        elif self.type=='log':
            if sc == 0:
                self.coefs = np.array([S(0,0)])

            else:
                N = 2 ** (sc - 1) + 1
                ag = np.linspace(-self.size, self.size, N)
                xg, yg = np.meshgrid(ag, ag)
                xg = np.reshape(xg, xg.size)
                yg = np.reshape(yg, yg.size)
                self.coefs = S(xg, yg)

            #print("LOG self.coefs=", self.coefs)
            self.coefs = np.log(self.coefs)
            #print("LOG self.coefs=", self.coefs)

        """
        elif self.type=='wav':
            nsc = self.DW.getnscales()
            N = 2**nsc
            ag = np.linspace(-self.size, self.size, N)
            xg, yg = np.meshgrid(ag, ag)
            vg = V(xg, yg)
            self.DW.setval(vg)
            ncoef = self.DW.gettotcoefscale(sc)
            m = np .empty(ncoef)
            for i in range(ncoef):
                m[i] = self.DW.getcoef(i)
        """

        self.coef2mod(sc, self.coefs)

        return self.coefs

    def coef2mod(self, sc, m):

        if self.type=='lin':
            if sc == 0: N = 2
            else:       N = 2 ** (sc - 1) + 1

            ag=np.linspace(-self.size, self.size, N)
            xg, yg = np.meshgrid(ag, ag)
            xg = np.reshape(xg, xg.size)
            yg = np.reshape(yg, yg.size)

            #print("coef2mod sc=",sc," N=",N," m=",len(m))
            if sc>0: vg = m
            else:    vg = np.full(N*N, fill_value=m)

        elif self.type=='log':
            if sc == 0: N = 2
            else:       N = 2 ** (sc - 1) + 1

            ag=np.linspace(-self.size, self.size, N)
            xg, yg = np.meshgrid(ag, ag)
            xg = np.reshape(xg, xg.size)
            yg = np.reshape(yg, yg.size)

            #print("coef2mod sc=",sc," N=",N," m=",len(m))
            if sc>0: vg = np.exp(m)
            else:    vg = np.full(N*N, fill_value=np.exp(m))

        """
        elif self.type=='wav':
            N = self.DW.base.shape[0]
            ag = np.linspace(-self.size, self.size, N)
            xg, yg = np.meshgrid(ag, ag)
            xg = np.reshape(xg, xg.size)
            yg = np.reshape(yg, yg.size)

            k=0
            for i in range(sc+1):
                ncoef = self.DW.getncoefscale(i)
                for j in range(ncoef):
                    self.DW.setcoefsc(i,j,m[k])
                    k=k+1

            vg = self.DW.getidwt()
            vg = np.reshape(vg, vg.size)
        """

        #print("coef2mod interp")
        #print(np.vstack((xg, yg)).T)
        #print(vg)
        self.S = LinearNDInterpolator(np.vstack((xg, yg)).T, vg)

        return self.S

    def check_valid(self, sc, m):

        if self.type=='lin':
            if sc == 0: N = 2
            else:       N = 2 ** (sc - 1) + 1

            ag=np.linspace(-self.size, self.size, N)
            xg, yg = np.meshgrid(ag, ag)
            xg = np.reshape(xg, xg.size)
            yg = np.reshape(yg, yg.size)

            #print("coef2mod sc=",sc," N=",N," m=",len(m))
            if sc>0: vg = m
            else:    vg = np.full(N*N, fill_value=m)

            if np.min(vg)<0.: return False

        return True

    def gettikhonov(self, sc):

        print("GETTIKHONOV type=",self.type)
        if self.type == 'lin' or self.type=='log':
            if sc == 0:
                N = 2
                L=np.array([1])
                return L

            else:
                N = 2 ** (sc - 1) + 1

            L = np.zeros((N**2,N**2))

            stencil=np.array([[-0.25,-0.5,-0.25], [-0.5,3,-0.5], [-0.25,-0.5,-0.25]])

            for i in range(N):
                for j in range(N):

                    r=i*N+j

                    for ii in range(3):
                        for jj in range(3):
                            ig=i+ii-1
                            jg=j+jj-1
                            if ig<0 or ig>N-1 or jg<0 or jg>N-1: continue
                            s=ig*N+jg
                            L[r,s]=stencil[ii,jj]

        #print("GETTIKH")
        #print(L)

        #L=L*(N/np.linalg.norm(L))

        return L

    def getnodes(self, sc):

        if self.type=='lin' or self.type=='log':
            if sc == 0: N = 2
            else:       N = 2 ** (sc - 1) + 1

            ag=np.linspace(-self.size, self.size, N)
            xg, yg = np.meshgrid(ag, ag)
            xg = np.reshape(xg, xg.size)
            yg = np.reshape(yg, yg.size)

        """
        elif self.type=='wav':
            N = self.DW.base.shape[0]
            ag = np.linspace(-self.size, self.size, N)
            xg, yg = np.meshgrid(ag, ag)
            xg = np.reshape(xg, xg.size)
            yg = np.reshape(yg, yg.size)

            k=0
            for i in range(sc+1):
                ncoef = self.DW.getncoefscale(i)
                for j in range(ncoef):
                    self.DW.setcoefsc(i,j,m[k])
                    k=k+1

            vg = self.DW.getidwt()
            vg = np.reshape(vg, vg.size)
        """

        return xg, yg

    def setInt(self, S):
        self.S=S

    def getV(self, x, y):
        return 1./self.S(x,y)
