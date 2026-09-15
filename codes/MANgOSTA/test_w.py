
import pywt
import numpy as np
import matplotlib.pyplot as plt

mat=np.zeros((32,32))
print(mat.shape)

wav='sym2'

cA, (cH, cV, cD)  = pywt.dwt2(mat, wav)
print(cA.shape, cH.shape, cV.shape, cD.shape)
mat=cA

cA, (cH, cV, cD)  = pywt.dwt2(mat, wav)
print(cA.shape, cH.shape, cV.shape, cD.shape)
mat=cA

cA, (cH, cV, cD)  = pywt.dwt2(mat, wav)
print(cA.shape, cH.shape, cV.shape, cD.shape)
mat=cA

cA, (cH, cV, cD)  = pywt.dwt2(mat, wav)
print(cA.shape, cH.shape, cV.shape, cD.shape)
mat=cA

cA, (cH, cV, cD)  = pywt.dwt2(mat, wav)
print(cA.shape, cH.shape, cV.shape, cD.shape)


"""
coef = pywt.wavedec2(mat, 'coif4')

print(len(coef))

coef[0][0][0]=1

m2=pywt.waverec2(coef, 'coif4')

plt.pcolor(m2)
plt.show()
"""