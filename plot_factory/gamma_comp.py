import math

import matplotlib.pyplot as plt


def gamma_paper(r):
    return 900/(r**8 + 3.584**8)

def gamma_fortran(r, p, c, s):
    return abs(-0.5*p*math.erf(s*(r - c)) - 1)

p = 0.0105
c = 2.8
s = 0.8

r = [*range(101)]

plt.plot(r, [gamma_paper(x) for x in r])
plt.plot(r, [gamma_fortran(x, p, c, s) for x in r])
plt.show()
plt.clf()
