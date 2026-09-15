
import numpy as np
from copy import deepcopy

class Inverse:

    def __init__(self, params):

        self.theta=params['regul']

    def settickhonov(self, L):
        #self.L = deepcopy(L)
        self.LL = np.dot(L.T,L)

    """
    def get_curv(self, p1, p2, p3):
        
        #Returns the center and radius of the circle passing the given 3 points.
        #In case the 3 points form a line, returns (None, infinity).
        
        temp = p2[0] * p2[0] + p2[1] * p2[1]
        bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
        cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
        det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])

        if abs(det) < 1.0e-6:
            print("EEE")
            return (None, np.inf)

        # Center of circle
        cx = (bc * (p2[1] - p3[1]) - cd * (p1[1] - p2[1])) / det
        cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det

        radius = np.sqrt((cx - p1[0]) ** 2 + (cy - p1[1]) ** 2)
        return 1. / radius
    """

    def getLcurve_sel(self, lm, ld):

        llm = np.log(lm)
        lld = np.log(ld)
        llm = (llm - np.min(llm)) / (np.max(llm) - np.min(llm))
        lld = (lld - np.min(lld)) / (np.max(lld) - np.min(lld))

        #kk = llm + lld
        #idx = np.argmin(kk)

        curv1 = np.zeros(len(llm))
        for i in range(1, len(llm) - 1):
            zeta1 = (lld[i + 1] - lld[i - 1]) / 2
            zeta2 = (lld[i + 1] + lld[i - 1] - 2 * lld[i])
            eta1 = (llm[i + 1] - llm[i - 1]) / 2
            eta2 = (llm[i + 1] + llm[i - 1] - 2 * llm[i])
            curv_num = zeta1 * eta2 - zeta2 * eta1
            #curv_den = np.power(zeta1 ** 2 + eta1 ** 2, 1.5)
            curv_den = (zeta1 ** 2 + eta1 ** 2) ** 1.5
            curv1[i] = -curv_num / curv_den

        # Check all the relative maxima
        rmax=[]
        cmax=int(len(curv1)*0.75)
        for i in range(1,cmax):
            if curv1[i]>curv1[i-1] and curv1[i]>curv1[i+1]:
                #print("Rel. maxim ",i)
                #if i<len(curv1)-3:
                rmax.append(i)

        if rmax==[]: idx=np.argmax(curv1)
        else:
            imax=0
            for i in range(len(rmax)):
                if curv1[rmax[i]]>curv1[rmax[imax]]:
                    imax=i
            idx = rmax[imax]

        #idx=np.argmax(curv1)

        #print("maximo ",idx, flush=True)

        """
        import matplotlib.pyplot as plt

        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(5, 1)

        ax = fig.add_subplot(gs[0, 0])
        ax.plot(curv1, 'k')
        for i in range(len(rmax)): ax.plot(rmax[i], curv1[rmax[i]], 'bo')
        ax.plot(idx, curv1[idx], 'ro')
        ax.set_title('curv')

        ax = fig.add_subplot(gs[1:, 0])
        ax.plot(lld, llm, 'ko')
        ax.plot(lld[0], llm[0], 'go')
        for i in range(len(rmax)): ax.plot(lld[rmax[i]], llm[rmax[i]], 'bo')
        ax.plot(lld[idx], llm[idx], 'ro')

        plt.show()
        """

        return idx

    def getinverse_tikh(self, G, m0, res, valid, inv, mref=None):

        print("Tikhonov regularization")

        M = G.shape[1]

        eps2r = np.flipud(np.logspace(-5, 10, 100))

        G1=np.dot(G.T,G)

        ####################
        lm = []
        ld = []
        gcv = []
        leps = []
        for i in range(len(eps2r)):

            eps0 = eps2r[i] * (1.-self.theta)
            eps2 = eps2r[i] * self.theta

            G2 = G1 + eps0*np.eye(M) + eps2*self.LL
            G3 = np.linalg.inv(G2)

            if inv=='LINEAR':
                G4 = np.dot(G3,G.T)
                dm = np.dot(G4, res)

            else:
                G4 = np.dot(G.T,res) + eps0 * (m0-mref) + eps2 * np.dot(self.LL, m0-mref)
                dm = np.dot(G3, G4)

            sy = np.dot(G, dm)
            m1 = m0 - dm
            if not valid(m1):
                #print(eps2r[i],"not valid ")
                break

            if inv=='LINEAR':
                lm.append(np.linalg.norm(dm))
            else:
                lm.append(np.linalg.norm(m1-mref))
            ld.append(np.linalg.norm(res-sy))

            gcv_num = np.linalg.norm(res-sy)**2
            Gig = np.dot(G3,G.T)
            A = np.dot(G,Gig)
            NA = A.shape[0]
            gcv_den = np.linalg.norm(np.diag(np.eye(NA)-A))**2
            gcv.append(gcv_num/gcv_den)

            leps.append(eps2r[i])

        gcv = np.array(gcv)
        idx1 = np.argmin(gcv)
        idx2 = self.getLcurve_sel(lm, ld)

        idx=min(idx1,idx2)
        eps2 = eps2r[idx]
        print("idx gcv=",idx1," idx L-curve=",idx2)
        print("best idx=%d eps2=%g" % (idx,eps2))

        eps0 = eps2 * (1. - self.theta)
        eps2 = eps2 * self.theta
        G2 = G1 + eps0 * np.eye(M) + eps2 * self.LL
        G3 = np.linalg.inv(G2)

        if inv == 'LINEAR':
            G4 = np.dot(G3, G.T)
            dm = np.dot(G4, res)

        else:
            G4 = np.dot(G.T, res) + eps0 * (m0 - mref) + eps2 * np.dot(self.LL, m0 - mref)
            dm = np.dot(G3, G4)

        return -dm, True

    def getinverse_tsvd(self, G, m0, res, valid, inv, mref=None):

        print("TSVD regularization")

        M = G.shape[1]

        # Truncated SVD
        U, S, VT = np.linalg.svd(G)
        V = VT.T

        maxp = len(S)

        if maxp > 1:

            lm = []
            ld = []
            lp = []
            #gcv = []
            for p in range(1, maxp + 1):

                Sp = S[range(0, p)]
                if Sp[-1] / Sp[0] < 1e-10: break

                Up = U[:, range(0, p)]
                Vp = V[:, range(0, p)]
                iSp = np.diag(1. / Sp)
                Gig = np.dot(np.dot(Vp, iSp), Up.T)
                dm = np.dot(Gig, res)
                sy = np.dot(G, dm)

                m1 = m0 - dm
                if not valid(m1): break

                if inv=='LINEAR':
                    lm.append(np.linalg.norm(dm))
                else:
                    lm.append(np.linalg.norm(m1-mref))

                ld.append(np.linalg.norm(res-sy))
                lp.append(p)

            if lp == []:
                return None, False

            idx = self.getLcurve_sel(lm, ld)
            popt = lp[idx]

        else:
            popt = 1

        print("popt=", popt)
        Up = U[:, range(0, popt)]
        Vp = V[:, range(0, popt)]
        Sp = S[range(0, popt)]
        iSp = np.diag(1. / Sp)
        Gig = np.dot(np.dot(Vp, iSp), Up.T)

        dm = np.dot(Gig, res)

        return -dm, True

    def getinverse(self, G, m0, res, valid, inv, mref=None):

        if self.theta==-1:
            dm, flag = self.getinverse_tsvd(G, m0, res, valid, inv, mref)
        else:
            dm, flag = self.getinverse_tikh(G, m0, res, valid, inv, mref)

        return dm, flag
