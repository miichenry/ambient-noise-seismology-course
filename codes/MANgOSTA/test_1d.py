
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def getarraycoord(Xp, Yp):

    def getpt(Xp, Yp, b, c):

        Xi = (b * (b * Xp - Yp) - c) / (1.+b**2)
        Yi = (-b * Xp + Yp - b * c)  / (1.+b**2)

        return Xi, Yi

    def adist(Xp, Yp, b, c):

        td = 0.
        for i in range(len(Xp)):
            Xi, Yi = getpt(Xp[i], Yp[i], b, c)
            td = td + np.sqrt( (Xi-Xp[i])**2 + (Yi-Yp[i])**2 )
            #td = td + (Xi-Xp[i])**2 + (Yi-Yp[i])**2

        return td

    def getbest_bc(Xp, Yp):

        mis = lambda m: adist(Xp, Yp, m[0], m[1])

        # Farthest couple
        Npts = len(Xp)
        dmax = 0.
        imax = 0
        jmax = 0
        for i in range(Npts):
            for j in range(i + 1):
                dist = np.sqrt((Xp[i] - Xp[j]) ** 2 + (Yp[i] - Yp[j]) ** 2)
                if dist > dmax:
                    dmax = dist
                    imax = i
                    jmax = j

        # Initial trial
        btrial = (Xp[jmax]-Xp[imax])/(Yp[imax]-Yp[jmax])
        ctrial = -(Xp[imax]+btrial*Yp[imax])

        m0 = np.array([btrial,ctrial])
        res = minimize(mis, m0, method='Nelder-Mead')

        #print(m0, res.x)

        b = res.x[0]
        c = res.x[1]

        return b, c

    ###############################

    Xc = np.average(Xp)
    Yc = np.average(Yp)

    Xpr = Xp - Xc
    Ypr = Yp - Yc

    b, c = getbest_bc(Xpr, Ypr)

    Npts = len(Xpr)
    Xi = np.zeros(Npts)
    Yi = np.zeros(Npts)

    for i in range(len(Xi)):
        Xi[i], Yi[i] = getpt(Xpr[i], Ypr[i], b, c)

    #print(Xi)
    #print(Yi)

    # Farthest couple
    dmax = 0.
    imax = 0
    jmax = 0
    for i in range(Npts):
        for j in range(i+1):
            dist = np.sqrt( (Xi[i]-Xi[j])**2 + (Yi[i]-Yi[j])**2 )
            if dist>dmax:
                dmax=dist
                imax=i
                jmax=j

    #print(imax, jmax)

    pos = np.zeros(Npts)
    for i in range(Npts):
        pos[i] = np.sqrt( (Xi[i]-Xi[imax])**2 + (Yi[i]-Yi[imax])**2 )

    return (Xi[imax]+Xc,Yi[imax]+Yc), (Xi[jmax]+Xc,Yi[jmax]+Yc), pos

#########################################
Xp=np.array([0, 2, 3, 5, 7, 10, 13])
Yp=np.array([-2, -1, -1, 1, 5, 4, 7])

P1, P2, pos = getarraycoord(Xp, Yp)

print(P1, P2)

plt.plot([P1[0],P2[0]], [P1[1],P2[1]], 'ko-')

plt.plot(Xp, Yp, 'rx')

plt.show()
