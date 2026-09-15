
import numpy as np

curv=np.zeros(100)
curv[30]=1.
sm = np.array([0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25])
curv_sm = np.convolve(curv, sm, mode='same')
print(np.argmax(curv), np.argmax(curv_sm))
exit(1)

"""
a=np.logspace(-1,1,10)
print(a)
b=np.flipud(a)
print(b)
"""

"""
from Tomo.velmod import Velmod

box=[0,0,10]
VM=Velmod(box)

sc=4

print("scale=",sc)
print("ncoef=",VM.getncoefsc(sc))
L=VM.gettikhonov(sc)
#print(L)

import matplotlib.pyplot as plt

plt.pcolor(L)
plt.show()
"""