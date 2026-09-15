
import numpy as np
import matplotlib.pyplot as plt

x,y=np.meshgrid(np.linspace(0.5,2.5,100),np.linspace(0.5,2.5,100))

k=np.sqrt(2)/2

ix=k*(x-y)
iy=k*(x+y)

iix=int((ix-0.5)/0.5)
iiy=int((iy-0.5)/0.5)
s=(iix+iiy) % 2

plt.pcolor(s)
plt.show()
