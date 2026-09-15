
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

N=5
x, y = np.meshgrid(np.linspace(-10,10,N), np.linspace(-10,10,N))
#x = np.linspace(-10,10,N)
#y = np.linspace(-10,10,N)
v=np.zeros((N,N))
v[2,2]=1.
#isp = interpolate.RectBivariateSpline(x,y,v)
x = np.reshape(x, x.size)
y = np.reshape(y, y.size)
v = np.reshape(v, v.size)
isp = interpolate.CloughTocher2DInterpolator(np.vstack((x, y)).T, v)
#isp = interpolate.LinearNDInterpolator(np.vstack((x, y)).T, v)

xnew, ynew = np.mgrid[-10:10:30j, -10:10:30j]
znew = isp(xnew, ynew)

plt.pcolor(znew)
plt.colorbar()
plt.show()
