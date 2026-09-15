
import numpy as np

def rot_compl(V):

    V=V/np.linalg.norm(V)

    alpha = np.arange(0, 180, 1) * np.pi / 180

    imax = 0
    Xmax = 0.
    for i in range(len(alpha)):
        cis = np.cos(alpha[i]) + np.sin(alpha[i]) * 1j
        Xx = np.real(V[0] * cis) ** 2
        Xy = np.real(V[1] * cis) ** 2
        Xz = np.real(V[2] * cis) ** 2
        X = np.sqrt(Xx + Xy + Xz)
        if X > Xmax:
            Xmax = X
            imax = i

    cis = np.cos(alpha[imax]) + np.sin(alpha[imax]) * 1j
    Xx = V[0] * cis
    Xy = V[1] * cis
    Xz = V[2] * cis

    U = np.array([Xx, Xy, Xz])
    U = U / np.linalg.norm(U)
    return np.real(U), Xmax

v=np.array([1+2j, -2+3j, 0-1j])
print(rot_compl(v))
