
import sys
import numpy as np
from Parser import Parser
from Tomo.tomo import Tomo
from version import getver
import gc

if __name__ == '__main__':

    print("MANGOSTA tomo ",getver())
    gc.disable()

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

    # ttcalc: STRAIGHT, SP
    params['ttcalc'] = 'STRAIGHT'

    params['ttgrid'] = 20

    # inversion: LINEAR, NONLINEAR
    params['inversion'] = 'LINEAR'

    # regul: 0 Tikhonov 0-th, 1 Tikhonov 2-th, -1 TSVD
    params['regul'] = 0.

    params['invdm'] = 1
    params['invmaxiter'] = 10

    while not P.finished():

        row = P.get()
        print(row)

        if row[0].upper() == "DISPDIR":
            params['dispdir'] = row[1]

        if row[0].upper() == "TOMODIR":
            params['tomodir'] = row[1]

        if row[0].upper() == "STAFILE":
            params['stafile'] = row[1]

        if row[0].upper() == "COMP":
            comp=[]
            for i in range(1,len(row)):
                comp.append(row[i])
            params['comp'] = comp

        if row[0].upper() == "TTCALC":
            ttcalc = row[1].upper()
            params['ttcalc'] = ttcalc
            if ttcalc=="SP": params['ttgrid'] = int(row[2])

        if row[0].upper() == "PERIODS":
            periods = []
            for i in range(len(row)-1):
                periods.append(float(row[i+1]))
            params['periods'] = periods

        if row[0].upper() == "BOX":
            params['box'] = row[1:]

        if row[0].upper() == "SCALE":
            params['scale'] = int(row[1])

        if row[0].upper() == "MULTISCALE":
            params['multiscale'] = int(row[1])

        if row[0].upper() == "INVERSION":
            params['inversion'] =  row[1].upper()
            count=2
            while count<len(row):
                if row[count].upper() == "DM":
                    params['invdm'] = float(row[count+1])
                    count = count + 2
                elif row[count].upper() == "MAXITER":
                    params['invmaxiter'] = int(row[count+1])
                    count = count + 2
                else:
                    print("Wrong input parameter ",row[count].upper())
                    exit(1)

        if row[0].upper() == "REGUL":
            if row[1].upper() == "TSVD":
                params['regul'] = -1.
            else:
                params['regul'] = float(row[1])

        if row[0].upper() == "DOTOMO":

            if not 'box' in params.keys():
                print("Error: parameter BOX required")
                exit(1)

            if not 'periods' in params.keys():
                print("Error: parameter PERIODS required")
                exit(1)

            if params['inversion']=="LINEAR" and params['ttcalc']=="SP":
                print("Error: Shortest path computation incompatible with linear inversion")
                exit(1)

            if params['inversion']=="LINEAR" and 'multiscale' in params.keys():
                print("Error: Multiscale incompatible with linear inversion")
                exit(1)

            for i in range(len(params['comp'])):

                print("COMP=",params['comp'])
                T=Tomo(params)
                T.dotomo()

        if row[0].upper() == "EXIT":
            exit(0)
