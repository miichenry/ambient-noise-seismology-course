# coding: utf-8

from Finder.finder import Finder
from Correlator.correlator import Correlator
from Censor.censor import Censor

print("Finder")
f=Finder('test', 'stacoord.dat', 'mseed')
print(f.stations)

f.check_windows(wpwr=16)

# Censor
print("Censor")
#ce=Censor(f)
#ce.check_netcoh()

print("Correlator")
c=Correlator(f, wav='sym8')
c.correlate(0.01,1.0)

"""
print "REPORT"
f.reportwin('out/win.log')
print "DONE"
"""

"""
nwin=len(f.twin)
for i in range(nwin):
    print f.twin[i]
    out=f.getallz(i)
    fname="out/%05d.npy" % i
    numpy.save(fname,out)
"""


"""
nw = f.getnwindows('TNOR', 'TSIS')
np = f.nptw
cc = numpy.zeros((nw,np))

for k in range(nw):
    w = f.getwindow('TNOR', 'TSIS', k)
    s1 = w[0] - numpy.mean(w[0])
    s2 = w[1] - numpy.mean(w[1])
    #cc[k,:] = cross_correlation_using_fft(s1, s2)
    #dw1=DWTO(s1)
    #dw1.getsigsc(2)
"""
