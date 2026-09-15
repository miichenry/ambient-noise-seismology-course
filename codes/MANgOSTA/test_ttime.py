
from Tomo.velmod import Velmod
from Tomo.ttime2D import TTime2D

T=TTime2D(20,20,1000)

V=Velmod([10,10,10000])
V.homo2coef(1,1000.)

T.setvelmod(V)

T.calc(5000,5000)
