# plots the multi-atom gamma

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf as erf

# from plot template
from matplotlib import rcParams
from matplotlib import ticker as mticker
from matplotlib.ticker import MultipleLocator
rcParams.update({'font.size': 32}) # changed from 36 to 32
rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})
rcParams['font.sans-serif'] = "Helvetica"
rcParams['font.family'] = "sans-serif"
rcParams['lines.linewidth'] = 2
    # width and length of the ticks
twdith = 2.5
rcParams.update({'xtick.major.width': twdith})
rcParams['ytick.major.width'] = twdith
rcParams['xtick.minor.width'] = twdith
rcParams['ytick.minor.width'] = twdith
tlenmaj = 9
rcParams['xtick.major.size'] = tlenmaj
rcParams['ytick.major.size'] = tlenmaj
tlenmaj = 6
rcParams['xtick.minor.size'] = tlenmaj
rcParams['ytick.minor.size'] = tlenmaj
rcParams['xtick.major.pad']='8'
rcParams['ytick.major.pad']='8'
    # colors
t4imaroon = '#941751'
t4iblue = '#006699'
t4ipetrol = '#007E71'
t4ipurp = '#531B93'
t4iorange = '#FF9300'
t4irose = '#FF818E'
t4ilime = '#8EFA00'
t4iturq = '#00FDFF'
t4istraw = '#FF2F92'
t4ivio = '#9E7BFF'
    # setup
fig = plt.figure(figsize=(10.5, 7))
ax = fig.add_subplot(111)
plt.setp(ax.spines.values(), lw=3)
fig.subplots_adjust(top=0.985, bottom=0.15, left=0.20, right=0.97, hspace=0.05, wspace=0.2) # changed from 0.17 to 20
#


def gamma_paper(r):
    return 900/(r**8 + 3.584**8)

def gamma_fortran(r, p, c, s):
    return np.abs(-0.5 * p * (erf(s*(r - c)) - 1))

# settings
folder = "C:/Users/essle/Downloads"
mode = "s"
p = 0.0205 #a
c = 2.8 #b
s = 0.8 #c
r = np.linspace(0.1, 10, 1000)
#

# plot
#plt.plot(r, gamma_paper(r), color=t4iblue, label =r"$\gamma^{theory}$")

if mode == "p":
    plt.plot(r, gamma_fortran(r, p, c, s), color=t4iblue, label =r"$\gamma_{p=0.0205, c=2.8, s=0.8}$")
    plt.plot(r, gamma_fortran(r, 0.0105, c, s), color=t4iorange, label =r"$\gamma_{p=0.0105, c=2.8, s=0.8}$")
    plt.plot(r, gamma_fortran(r, 0.0005, c, s), color=t4imaroon, label =r"$\gamma_{p=0.0005, c=2.8, s=0.8}$")
elif mode == "c":
    plt.plot(r, gamma_fortran(r, 0.0105, 3.8, s), color=t4iblue, label =r"$\gamma_{p=0.0105, c=3.8, s=0.8}$")
    plt.plot(r, gamma_fortran(r, 0.0105, 2.8, s), color=t4iorange, label =r"$\gamma_{p=0.0105, c=2.8, s=0.8}$")
    plt.plot(r, gamma_fortran(r, 0.0105, 1.8, s), color=t4imaroon, label =r"$\gamma_{p=0.0105, c=1.8, s=0.8}$")
elif mode == "s":
    plt.plot(r, gamma_fortran(r, 0.0105, 2.8, 1.3), color=t4iblue, label =r"$\gamma_{p=0.0105, c=2.8, s=1.3}$")
    plt.plot(r, gamma_fortran(r, 0.0105, 2.8, 0.8), color=t4iorange, label =r"$\gamma_{p=0.0105, c=2.8, s=0.8}$")
    plt.plot(r, gamma_fortran(r, 0.0105, 2.8, 0.3), color=t4imaroon, label =r"$\gamma_{p=0.0105, c=2.8, s=0.3}$")

plt.xlabel(r"r")
plt.ylabel(r"$\gamma$")

plt.legend(fontsize=20, frameon=False)
plt.savefig(os.path.join(folder, f"gamma_comp_{mode}.pdf"), dpi=600)
plt.clf()

    # from plot template
ax.xaxis.set_major_locator(MultipleLocator(10)) # changed from 5 to 10
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
#
