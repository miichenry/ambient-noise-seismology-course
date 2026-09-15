
import sys, os
import numpy as np
import utm
from scipy.interpolate import LinearNDInterpolator
from Parser import Parser
from Station.station import Station
from version import getver

from Tomo.ttime2D import TTime2D
from Tomo.velmod import Velmod

if __name__ == '__main__':

    print("MANGOSTA synth ",getver())

    if len(sys.argv)<2:
        raise NameError('No input file specified!')

    try:
        nlines = sum(1 for line in open(sys.argv[1]))
    except:
        raise NameError('Wrong input file specified!')

    P=Parser(sys.argv[1])

    ###############################
    # Interpretation of input file

    # Default parameters
    params = {}

    params['dispdir'] = 'disp'
    params['tomodir'] = 'tomo'
    params['stafile'] = 'stations.csv'
    params['comp'] = ['ZZ']
    params['scale'] = 1

    params['ttcalc'] = 'STRAIGHT'
    params['ttgrid'] = 20

    params['sigma'] = 0.0

    while not P.finished():

        row = P.get()
        print(row)

        if row[0].upper() == "DISPDIR":
            params['dispdir'] = row[1]

        if row[0].upper() == "STAFILE":
            params['stafile'] = row[1]

        if row[0].upper() == "COMP":
            comp = []
            for i in range(1, len(row)):
                comp.append(row[i])
            params['comp'] = comp

        if row[0].upper() == "PERIODS":
            periods = []
            for i in range(len(row)-1):
                periods.append(float(row[i+1]))
            params['periods'] = periods

        if row[0].upper() == "BOX":
            params['box'] = row[1:]

        if row[0].upper() == "MODEL":
            params['model'] = row[1]

        if row[0].upper() == "TTCALC":
            ttcalc = row[1].upper()
            params['ttcalc'] = ttcalc
            if ttcalc == "SP": params['ttgrid'] = int(row[2])

        if row[0].upper() == "SIGMA":
            params['sigma'] = float(row[1])

    if not 'box' in params.keys():
        print("Error: parameter BOX required")
        exit(1)

    if not 'periods' in params.keys():
        print("Error: parameter PERIODS required")
        exit(1)

    if not 'model' in params.keys():
        print("Error: parameter MODEL required")
        exit(1)

    st = Station(params['stafile'])

    # Compute stations coordinates in box reference system
    nsta=len(st.stalat)
    LAS=st.stalat
    LOS=st.stalon
    STNA=st.staname

    lat0 = float(params['box'][0])
    lon0 = float(params['box'][1])
    size = float(params['box'][2])
    X0, Y0, Z0, L0 = utm.from_latlon(lat0, lon0)

    XS = np.empty(nsta)
    YS = np.empty(nsta)
    for i in range(nsta):
        xs, ys, zs, ls = utm.from_latlon(LAS[i], LOS[i], Z0)
        XS[i] = xs - X0
        YS[i] = ys - Y0
        if XS[i]<-size or XS[i]>size or YS[i]<-size or YS[i]>size:
            print("Error: station ",STNA[i]," outside box")
            exit(1)

    # Read model
    mod = np.loadtxt(params['model'])
    LAP = mod[:, 0]
    LOP = mod[:, 1]
    VP = mod[:, 2] * 1000.
    print("VP ",np.min(VP),np.max(VP))

    # Interpolate model
    XP = np.empty(LAP.shape)
    YP = np.empty(LAP.shape)
    for i in range(len(XP)):
        XP[i], YP[i], zs, ls = utm.from_latlon(LAP[i], LOP[i], Z0)

    XP = XP - X0
    YP = YP - Y0
    C=np.vstack((XP, YP)).T
    V = LinearNDInterpolator(C, VP)

    if np.min(XP)>-size or np.max(XP)<size or np.min(YP)>-size or np.max(YP)<size:
        print("Error: input model does not cover the box")
        exit(1)

    print("X ", np.min(XP), np.max(XP))
    print("Y ", np.min(YP), np.max(YP))

    # Create dispdir
    if not os.path.isdir(params['dispdir']):
        os.mkdir(params['dispdir'])

    if params['ttcalc']=='STRAIGHT':

        NINT = 100
        for cmp in params['comp']:
            # Loop over pairs
            for i in range(nsta):
                if XS[i]<-size or XS[i]>size or YS[i]<-size or YS[i]>size: continue
                for j in range(i+1,nsta):
                    if XS[j] < -size or XS[j] > size or YS[j] < -size or YS[j] > size: continue
                    dist = np.sqrt( (XS[i]-XS[j])**2 + (YS[i]-YS[j])**2 )

                    XG = np.linspace(XS[i], XS[j], NINT)
                    YG = np.linspace(YS[i], YS[j], NINT)
                    VG = V(XG, YG)
                    DS = dist/NINT
                    DELTAT = DS / VG
                    TG = np.sum(DELTAT)

                    fname = "%s/disp_%s_%s_%s_sy.dat" % (params['dispdir'],cmp,STNA[i],STNA[j])
                    f=open(fname,'w')
                    f.write("%f\n" % (dist/1000.) ) #!! Dist in km
                    for p in params['periods']:
                        vel = 0.001*dist/TG + np.random.randn() * params['sigma']
                        f.write("%f %f -1\n" % (p,vel)) #!! VEL in km/s
                    f.close()

    else:

        # Create velocity model and TTgrid
        DL = (size * 2) / params['ttgrid']
        xl = (np.arange(params['ttgrid']) + 0.5) * DL - size
        xg, yg = np.meshgrid(xl, xl)
        VM=Velmod(params['box'])
        V = LinearNDInterpolator(C, 1./VP)
        VM.setInt(V)

        print("ttgrid=",params['ttgrid']," DL=",DL," size=",size)
        TT2D = TTime2D(params['ttgrid'], params['ttgrid'], DL, -size, -size)
        TT2D.setvelmod(VM)

        for cmp in params['comp']:
            # Loop over sta1
            for i in range(nsta):
                if XS[i]<-size or XS[i]>size or YS[i]<-size or YS[i]>size: continue

                print(i,STNA[i])
                TT2D.calc(XS[i], YS[i])

                # Loop over sta2
                for j in range(i+1,nsta):
                    if XS[j] < -size or XS[j] > size or YS[j] < -size or YS[j] > size: continue
                    dist = np.sqrt((XS[i] - XS[j]) ** 2 + (YS[i] - YS[j]) ** 2)

                    TG=TT2D.gettime(XS[j], YS[j])
                    #print(j,TG)
                    fname = "%s/disp_%s_%s_%s_sy.dat" % (params['dispdir'], cmp, STNA[i], STNA[j])
                    f = open(fname, 'w')
                    f.write("%f\n" % (dist / 1000.))  # !! Dist in km
                    for p in params['periods']:
                        vel = 0.001 * dist / TG + np.random.randn() * params['sigma']
                        f.write("%f %f -1\n" % (p, vel))  # !! VEL in km/s
                    f.close()
