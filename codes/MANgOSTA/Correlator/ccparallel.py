#!/usr/bin/env python

import time
from mpi4py import MPI
import os
from Correlator.correlator import Correlator

def parallel(sta, f, ccdir, figdir, wavelet, fmin, fmax, cmp, smooth_wh, onebit):

    c = Correlator(f, ccdir, figdir, wav=wavelet)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    FLAG_AVAIL=0
    FLAG_ON=1
    FLAG_OFF=-1

    TAG_AVAIL=11
    TAG_FLAG=12
    TAG_DATA=13

    irecvsource = False

    # MASTER
    if rank == 0:

        count = 0

        # Request staus flag
        req = [True] * (size - 1)

        # INIT

        for k in range(1, size):
            if irecvsource:
                req[k - 1] = comm.irecv(source=k, tag=TAG_AVAIL)  # recibe una solicitud
            else:
                req[k - 1] = comm.irecv(dest=k, tag=TAG_AVAIL)  # recibe una solicitud

        start = time.time()

        while count < P.num_points:

            for k in range(1, size):

                if req[k - 1].Test():
                    # print "%s process %d available" % (time.strftime('%H:%M:%S'), k)
                    rflag = req[k - 1].wait()

                    # print "%s process %d sending count=%d x,y,z=%f,%f,%f lat,lon,dep=%f,%f,%f" % (time.strftime('%H:%M:%S'), k, count, x, y, z, lat, lon, dep)
                    comm.send(FLAG_ON, dest=k, tag=TAG_FLAG)
                    comm.send(sta, dest=k, tag=TAG_DATA)

                    count = count + 1

                    #if count % 1000 == 0:
                    #    elaps = time.time() - start
                    #    remain = elaps * float(P.num_points) / count - elaps

                    #if count >= P.num_points: break

                    if irecvsource:
                        req[k - 1] = comm.irecv(source=k, tag=TAG_AVAIL)
                    else:
                        req[k - 1] = comm.irecv(dest=k, tag=TAG_AVAIL)

                # else:
                #    print "%s process %d unavailable" % (time.strftime('%H:%M:%S'),k)

            # print "%s PAUSE\n" % (time.strftime('%H:%M:%S'))
            time.sleep(0.001)

        # Wait for pending slave processes
        flag_pending = True
        while flag_pending:

            flag_pending = False
            for k in range(1, size):
                if not req[k - 1].Test(): flag_pending = True
            time.sleep(1)

        # Terminate slaves
        for k in range(1, size):
            req[k - 1].wait()
            comm.send(FLAG_OFF, dest=k, tag=TAG_FLAG)

        MPI.Finalize()

    # SLAVE
    else:
        ful = open('out/ful_%d.out' % rank, 'w')
        sum = open('out/sum_%d.out' % rank, 'w')
        log = open('out/log_%d.out' % rank, 'w')

        log.write("rank %d starting init at %s\n" % (rank,time.strftime('%H:%M:%S')))

        c.correlate_trad(fmin, fmax, components=cmp, smooth_wh=smooth_wh, onebit_p=onebit)  #insertar la función correlate

        log.write("rank %d terminated init at %s\n" % (rank,time.strftime('%H:%M:%S')))

        while True:
             #log.write("rank %d waiting for data at %s\n" % (rank,time.strftime('%H:%M:%S')) )
            sreq = comm.isend(FLAG_AVAIL, dest=0, tag=TAG_AVAIL)

            sreq.wait()

            flag = comm.recv(source=0, tag=TAG_FLAG)

            if flag==FLAG_ON:
                #log.write('receiving data at %s\n' % time.strftime('%H:%M:%S'))
                x = comm.recv(source=0, tag=TAG_DATA)

            else:
                log.write("rank %d received KILL at %s\n" % (rank, time.strftime('%H:%M:%S')))
                log.close()
                ful.close()
                sum.close()
                MPI.Finalize()
                break