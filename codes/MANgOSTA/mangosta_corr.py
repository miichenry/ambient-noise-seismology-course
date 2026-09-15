#!/usr/bin/python3

# coding: utf-8

import sys, pickle

from Parser import Parser
from Station.station import Station
from Finder.finder import Finder
from Censor.censor import Censor
from Correlator.correlator import Correlator
#from FTAN.ftan import FTAN
from version import getver

if __name__ == '__main__':

    print("MANGOSTA correlator ", getver())

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

    wavelet = 'sym8'

    fmin = 0.1
    fmax = 1.0

    wlen = 200  # in seconds
    wpwr = -1

    datadir = 'data'
    dataext = ''
    stafile = 'stations.csv'
    st = []

    # Finder
    f = []

    ccdir = 'cc'
    figdir = 'fig'
    smooth_wh = 1   # 1 means no whitening
    onebit = 1.     # 1 means no onebit

    # FLAGS
    flagFinder = False
    flagCensor = False
    flagCorrelator = False
    flagFTAN = False

    print("START PROCESSING")

    while not P.finished():

        row = P.get()
        print(row)

        if row[0].upper() == "WAVELET":
            wavelet = row[1]

        elif row[0].upper() == 'FMIN':
            fmin = float(row[1])

        elif row[0].upper() == 'FMAX':
            fmax = float(row[1])

        elif row[0].upper() == 'WPWR':
            wpwr = int(row[1])

        elif row[0].upper() == 'WLEN':
            wlen = int(row[1])

        elif row[0].upper() == 'WHITENING':
            if row[1].upper() == 'FALSE':
                smooth_wh = 1
            elif row[1].upper() == 'UNIFORM':
                smooth_wh = -1
            else:
                smooth_wh = int(row[1])

        elif row[0].upper() == 'ONEBIT':
            if row[1].upper() == 'TRUE':
                onebit = 0
            if row[1].upper() == 'FALSE':
                onebit = 1
            else:
                onebit = float(row[1])

        elif row[0].upper() == 'DATADIR':
            datadir = row[1]

        elif row[0].upper() == 'DATAEXT':
            dataext = row[1]

        elif row[0].upper() == 'STAFILE':
            stafile = row[1]

        elif row[0].upper() == 'CCDIR':
            ccdir = row[1]

        elif row[0].upper() == 'FIGDIR':
            figdir = row[1]

        elif row[0].upper() == 'SAVE':

            if row[1].upper() == 'FINDER':

                if not flagFinder:
                    raise NameError('ERROR: SAVE before defining Finder')

                fh = open(row[2], 'wb')
                pickle.dump(f, fh, pickle.HIGHEST_PROTOCOL)
                fh.close()

        elif row[0].upper() == 'LOAD':

            if row[1].upper() == 'FINDER':

                fh = open(row[2], 'rb')
                f = pickle.load(fh)
                fh.close()

                flagFinder = True

        elif row[0].upper() == 'FINDER':
            print("Finder")

            st = Station(stafile)

            #f=Finder(datadir, st, dataext)
            f = Finder(datadir, st)
            f.check_windows(wlen_=wlen, wpwr_=wpwr)
            flagFinder = True
            #
            f.reportwin('win.log')

        elif row[0].upper() == 'NETCOH':

            if not flagFinder:
                raise NameError('ERROR: Starting NetCoh without doing Finder')

            print("Doing NetCoh")

            thr = float(row[1])

            fstack_min = float(row[2])
            fstack_max = float(row[3])

            if len(row)>4:
                subw_len = int(row[4])
                subw_n = int(row[5])
            else:
                subw_len = 30
                subw_n = 10

            #subnet = row[5].upper()

            cen = Censor(f, fmin, fmax)
            cen.check_netcoh(fstack_min, fstack_max, thr, subw_len, subw_n)

        elif row[0].upper() == 'CORRELATOR':
            
            if not flagFinder:
                raise NameError('ERROR: Starting Correlator without doing Finder')

            if len(row)<2:
                raise NameError('ERROR: please specify components in correlator')

            print("Correlator")
            c=Correlator(f, ccdir, figdir, wav=wavelet)
            #print("row1=",row[1].upper())
            if row[1].upper()=="ALL":
                cmp=['ZZ', 'ZR', 'ZT', 'RR', 'RT', 'TT']
            else:
                cmp=row[1:]

            #print("cmp=",cmp)
            c.correlate(fmin, fmax, components=cmp, smooth_wh=smooth_wh, onebit_p=onebit)

        elif row[0].upper() == 'CORRELATOR_TRAD':

            if not flagFinder:
                raise NameError('ERROR: Starting Correlator_Trad without doing Finder')

            if len(row)<2:
                raise NameError('ERROR: please specify components in correlator')

            print("Correlator")
            c=Correlator(f, ccdir, figdir, wav=wavelet)
            if row[1].upper()=="ALL":
                cmp=['ZZ', 'ZR', 'ZT', 'RR', 'RT', 'TT']
            else:
                cmp=row[1:]

            c.correlate_trad(fmin, fmax, components=cmp, smooth_wh=smooth_wh, onebit_p=onebit)

        elif row[0].upper() == 'FTAN':

            if not flagSTA:
                st = Station(stafile)
                flagSTA = True

            print("FTAN")
            ft=FTAN(ccdir, figdir, st)
            ft.doFTAN(fig=True)
            flagFTAN = True

        elif row[0].upper() == 'EXIT':
            print("END PROCESSING")
            exit(0)

        else:
            print("ERROR: wrong command ", row[0].upper())

    print("END PROCESSING")
