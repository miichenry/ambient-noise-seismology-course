
import numpy as np
from rdt import getrdt

lomin=2
lomax=3
lamin=1
lamax=2
nlon=10
nlat=10

NSTA=20
stalat=np.random.uniform(1.,2., NSTA)
stalon=np.random.uniform(2.,3., NSTA)

rdt, rd = getrdt(stalat, stalon, lomin, lomax, lamin, lamax, nlon, nlat, plot=True, scale=0.02)
