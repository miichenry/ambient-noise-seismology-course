
import numpy as np

a=np.eye(4)

a=np.array([[2,-1,0,0],
            [-1,2,-1,0],
            [0,-1,2,-1],
            [0,0,-1,2]])

a=a/np.linalg.norm(a)
print(np.linalg.norm(a)*2)
print(a)
