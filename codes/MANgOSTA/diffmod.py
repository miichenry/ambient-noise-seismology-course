
import sys
import numpy as np

def mad(v):
    #med=np.median(v)
    med = np.mean(v)
    res=np.abs(v-med)
    #print("median=", med, " mad=")
    return np.median(res)

data = sys.stdin.readlines()

mo=[]
mt=[]

for line in data:
    row=line.split(' ')
    mo.append(float(row[0]))
    mt.append(float(row[1]))

mo=np.array(mo)
mt=np.array(mt)

L2=np.linalg.norm(mo-mt,2)
L2_t=100*L2/np.linalg.norm(mt,2)

#L1=np.linalg.norm(mo-mt,1)
#L1_t=100*L1/np.linalg.norm(mt,1)

#Linf=np.max(np.abs(mo-mt))
#Linf_t=100*Linf/np.max(mt)
#MAD=mad(mo-mt)
#MAD_t=100*MAD/mad(mt)

#CC=np.sum(mo*mt)
#CC_t=100*CC/np.sum(mt*mt)

#print(L2,L1,MAD,L2_t,L1_t,MAD_t)
print(L2_t)
