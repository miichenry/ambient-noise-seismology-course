
from Finder.finder import Finder
from Censor.coherence import check_window_coherence

import psutil
import numpy as np

class Censor:

    f=[]    # Finder

    def __init__(self, f, fmin, fmax):

        if not isinstance(f,Finder):
            print("INTERNAL ERROR: init Censor without Finder!")
            raise Exception

        self.f = f
        self.fmin = fmin
        self.fmax = fmax

    def check_netcoh(self, fstack_min, fstack_max, coh_thr, subw_len, subw_n):

        swin=self.f.swin
        twin=self.f.twin

        nwin=swin.shape[0]
        nsta=swin.shape[1]

        wlen = self.f.nptw * self.f.delta

        log=open("netcoh.log","w")

        for i in range(nwin):

            strz=self.f.getStreamZ(i)

            test=check_window_coherence(starttime=twin[i], timedelta_hours=wlen/3600,
                                        ffilter_min=self.fmin, ffilter_max=self.fmax,
                                        fstack_min=fstack_min, fstack_max=fstack_max,
                                        subwindow_seconds=subw_len,
                                        subwindow_number=subw_n,
                                        work_path='',
                                        stream=strz,
                                        coherence_threshold=coh_thr,
                                        fig_coherencestack=False)

            print(twin[i]," ",test)
            log.write(twin[i].strftime("%Y-%m-%d %H:%M:%S")+" "+str(test)+"\n")
            log.flush()
            if test > coh_thr:
                print("Deleting window ",twin[i])
                log.write("Deleting window "+twin[i].strftime("%Y-%m-%d %H:%M:%S")+"\n")
                log.flush()
                for j in range(nsta):
                    swin[i,j]=False

            print(psutil.virtual_memory())

        log.close()
