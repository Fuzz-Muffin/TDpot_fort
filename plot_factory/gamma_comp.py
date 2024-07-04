import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf as erf


def gamma_paper(r):
    return 900/(r**8 + 3.584**8)

def gamma_fortran(r, p, c, s):
    return np.abs(-0.5 * p * (erf(s*(r - c)) - 1))

p = 0.0105 #a
c = 2.8 #b
s = 0.8 #c
r = np.linspace(0.1, 20, 1000)

plt.plot(r, gamma_paper(r))
plt.plot(r, gamma_fortran(r, p, c, s))
plt.show()
plt.clf()
