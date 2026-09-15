
import numpy as np
import matplotlib.pyplot as plt

def get_curv(p1, p2, p3):
    """
    Returns the center and radius of the circle passing the given 3 points.
    In case the 3 points form a line, returns (None, infinity).
    """
    temp = p2[0] * p2[0] + p2[1] * p2[1]
    bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
    cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
    det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])

    if abs(det) < 1.0e-6:
        print("EEE")
        return (None, np.inf)

    # Center of circle
    cx = (bc*(p2[1] - p3[1]) - cd*(p1[1] - p2[1])) / det
    cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det

    radius = np.sqrt((cx - p1[0])**2 + (cy - p1[1])**2)
    return 1./radius

"""
ang=np.linspace(0+0.1,np.pi/2-0.1,20)
a0=0
e=1.5
d=1
r=(e*d)/(1.-e*np.cos(ang-a0))
lld=r*np.sin(ang)
llm=r*np.cos(ang)
"""

lld=np.linspace(1.1,1.9,10)
#llm=(lld-0.5)**2
llm=np.sqrt(lld**2-1)

curv1 = np.zeros(len(llm))
for i in range(1, len(llm) - 1):
    zeta1 = (lld[i + 1] - lld[i - 1]) / 2
    zeta2 = (lld[i + 1] + lld[i - 1] - 2 * lld[i])
    eta1 = (llm[i + 1] - llm[i - 1]) / 2
    eta2 = (llm[i + 1] + llm[i - 1] - 2 * llm[i])
    curv_num = zeta1 * eta2 - zeta2 * eta1
    #curv_den = np.power(zeta1 ** 2 + eta1 ** 2, 1.5)
    curv_den = (zeta1 ** 2 + eta1 ** 2) ** 1.5
    curv1[i] = curv_num / curv_den

curv2 = np.zeros(len(llm))
for i in range(1, len(llm) - 1):
    h1 = lld[i] - lld[i - 1]
    h2 = lld[i+1] - lld[i]
    f0 = llm[i-1]
    f1 = llm[i]
    f2 = llm[i+1]
    y1 = -h2*f0/(h1*(h1+h2)) -(h1-h2)*f1/(h1*h2) + h1*f2/(h2*(h1+h2))
    y2 = 2*(h2*f0-(h1+h2)*f1+h1*f2)/(h1*h2*(h1+h2))
    curv2[i] = y2/np.power(1+y1**2, 1.5)

curv3 = np.zeros(len(llm))
for i in range(1, len(llm) - 1):
    p1 = np.array([lld[i - 1], llm[i - 1]])
    p2 = np.array([lld[i], llm[i]])
    p3 = np.array([lld[i + 1], llm[i + 1]])
    curv3[i]=get_curv(p1,p2,p3)

idx=np.argmax(curv1)

import matplotlib.pyplot as plt

fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(5, 1)

ax = fig.add_subplot(gs[0, 0])
ax.plot(curv1, 'k')
ax.plot(curv2, 'r')
ax.plot(curv3, 'g')
ax.plot(idx, curv1[idx], 'ro')
ax.set_title('curv')

ax = fig.add_subplot(gs[1, 0])
ax.plot(lld, 'k')
ax.plot(idx, lld[idx], 'ro')
ax.set_title('lld')

ax = fig.add_subplot(gs[2, 0])
ax.plot(llm, 'k')
ax.plot(idx, llm[idx], 'ro')
ax.set_title('llm')

ax = fig.add_subplot(gs[3:, 0])
ax.plot(lld, llm, 'ko')
ax.plot(lld[0], llm[0], 'bo')
ax.plot(lld[idx], llm[idx], 'ro')

plt.show()
