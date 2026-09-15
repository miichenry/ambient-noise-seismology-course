
import numpy as np
import utm
import matplotlib.pyplot as plt
from copy import deepcopy

# Intersection of line with rectangle
def liang_barsky_clipper(xmin,  ymin,  xmax,  ymax, x1,  y1,  x2,  y2):

    # defining variables
    p1 = -(x2 - x1)
    p2 = -p1
    p3 = -(y2 - y1)
    p4 = -p3

    q1 = x1 - xmin
    q2 = xmax - x1
    q3 = y1 - ymin
    q4 = ymax - y1

    posarr = np.zeros(5)
    negarr = np.zeros(5)

    posind = 1
    negind = 1
    posarr[0] = 1
    negarr[0] = 0

    if ((p1 == 0 and q1 < 0) or (p3 == 0 and q3 < 0)) :
        return []

    if (p1 != 0):
        r1 = q1 / p1
        r2 = q2 / p2
        if (p1 < 0):
            negarr[negind] = r1   # for negative p1, add it to negative array
            posarr[posind] = r2   # and add p2 to positive array
            negind = negind + 1
            posind = posind + 1
        else:
            negarr[negind] = r2
            posarr[posind] = r1
            negind = negind + 1
            posind = posind + 1

    if (p3 != 0):
        r3 = q3 / p3
        r4 = q4 / p4
        if (p3 < 0):
            negarr[negind] = r3
            posarr[posind] = r4
            negind = negind + 1
            posind = posind + 1
        else:
            negarr[negind] = r4
            posarr[posind] = r3
            negind = negind + 1
            posind = posind + 1

    rn1 = np.max(negarr[0:negind])  # maximum of negative array
    rn2 = np.min(posarr[0:posind])  # minimum of positive array

    if (rn1 > rn2):
        return None

    xn1 = x1 + p2 * rn1
    yn1 = y1 + p4 * rn1 # computing new points

    xn2 = x1 + p2 * rn2
    yn2 = y1 + p4 * rn2

    return [xn1, yn1, xn2, yn2]

def getrdt(stalat, stalon, lomin, lomax, lamin, lamax, nlon, nlat, plot=False, scale=0.01):

    # UTM mesh
    lat0 = (lamin + lamax) / 2
    lon0 = (lomin + lomax) / 2
    x0, y0, z0, l0 = utm.from_latlon(lat0, lon0)

    long, latg = np.meshgrid(np.linspace(lomin,lomax,nlon+1),np.linspace(lamin,lamax,nlat+1))

    xg = np.empty(long.shape)
    yg = np.empty(long.shape)
    for i in range(xg.shape[0]):
        for j in range(xg.shape[1]):
            xij, yij, zij, sij = utm.from_latlon(latg[i,j],long[i,j],z0)
            xg[i, j] = xij
            yg[i, j] = yij

    xg = xg - x0
    yg = yg - y0

    # Area
    ag = np.empty((nlat,nlon))

    for i in range(ag.shape[0]):
        for j in range(ag.shape[1]):
            ag[i,j]=(xg[i+1,j]-xg[i,j])*(yg[i,j+1]-yg[i,j])

    # UTM sta
    NSTA=stalat.size
    stax = np.empty(NSTA)
    stay = np.empty(NSTA)
    for i in range(NSTA):
        xs, ys, zs, ls = utm.from_latlon(stalat[i], stalon[i], z0)
        stax[i] = xs - x0
        stay[i] = ys - y0

    # Ray Density Tensor
    rdt=np.zeros((nlat,nlon,2,2), dtype=float)
    # Ray density (trace of RDT)
    rd = np.zeros((nlat, nlon), dtype=float)

    T0 = np.array([[1,0],[0,0]])

    # Loop over station pairs
    for i in range(stax.size):
        for j in range(i+1,stax.size):

            # Loop over grid cells
            for a in range(nlat):
                for b in range(nlon):
                    res = liang_barsky_clipper(xg[a,b],  yg[a,b],  xg[a,b+1],  yg[a+1,b], stax[i],  stay[i],  stax[j],  stay[j])
                    if res != None:
                        # Length
                        ls = np.sqrt( (res[2]-res[0])**2 + (res[3]-res[1])**2 )
                        # Azimuth
                        az = np.arctan2(res[3]-res[1],res[2]-res[0])
                        R=np.array([[np.cos(az),-np.sin(az)],[np.sin(az),np.cos(az)]])
                        # Tensor rotation and normalization
                        Tr=np.dot(R,np.dot(T0,R.T)) * ls/ag[a,b]
                        rdt[a,b,:,:] = rdt[a,b,:,:] + Tr
                        #rd[a,b]=1.

    for a in range(nlat):
        for b in range(nlon):
            rd[a,b]=rdt[a,b,0,0]+rdt[a,b,1,1]

    if plot==True:

        rdn = deepcopy(rd)
        rdn[rdn == 0] = np.nan

        plt.pcolor(long, latg, rd)
        for i in range(NSTA):
            for j in range(i, NSTA):
                plt.plot([stalon[i], stalon[j]], [stalat[i], stalat[j]], '-', color='gray')
        plt.plot(stalon, stalat, 'k*')

        # Plot RDT
        ang = np.linspace(-np.pi, np.pi, 100)
        xr = np.cos(ang)
        yr = np.sin(ang)
        pt = np.vstack((xr, yr))
        for a in range(nlat):
            for b in range(nlon):
                if np.isnan(rd[a, b]): continue
                Rab = np.squeeze(rdt[a, b, :, :])
                ptt = np.dot(Rab, pt)
                lon0 = (long[a, b] + long[a, b + 1]) / 2.
                lat0 = (latg[a, b] + latg[a + 1, b]) / 2.
                plt.plot(lon0 + ptt[0, :] * scale / rd[a, b], lat0 + ptt[1, :] * scale / rd[a, b], 'k')

        plt.show()

    return rdt, rd
