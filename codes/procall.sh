#!/bin/bash

DIR=/media/ivan/42922B2048EFB428/SCRIPTS/NOISETOMO/MANgOSTA_v2
#DIR=MANgOSTA


for regul in 0.50
do

for per in 0.9 1.0  1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0  2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0  4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0
do

### INVERT
ftomo="tomo.in"

echo "STAFILE stations.csv" > $ftomo
echo "COMP ZZ" >> $ftomo
echo "BOX  47.090751 8.272587 25000" >> $ftomo
echo "PERIODS 0.9 1.0  1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0  4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0" >> $ftomo
echo "PERIODS ${per}" >> $ftomo
echo "DISPDIR disp" >> $ftomo
echo "REGUL ${regul}" >> $ftomo
"""
echo "TTCALC STRAIGHT" >> $ftomo
echo "TOMODIR tomo_lin-$regul" >> $ftomo
echo "INVERSION LINEAR" >> $ftomo
echo "SCALE 7" >> $ftomo
echo "DOTOMO" >> $ftomo

echo "TOMODIR tomo_nolin-$regul" >> $ftomo
echo "TTCALC SP 100" >> $ftomo
echo "INVERSION NONLINEAR MAXITER 20 DM 0.1" >> $ftomo
echo "SCALE 7" >> $ftomo
echo "DOTOMO" >> $ftomo
"""
echo "TOMODIR tomo_multi-$regul" >> $ftomo
echo "TTCALC SP 100" >> $ftomo
echo "INVERSION NONLINEAR MAXITER 20 DM 1" >> $ftomo
echo "MULTISCALE 7" >> $ftomo
echo "DOTOMO" >> $ftomo

python3 ${DIR}/mangosta_tomo.py $ftomo

done

done

