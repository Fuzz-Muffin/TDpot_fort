import os

import matplotlib.pyplot as plt

# from plot template
from matplotlib import rcParams
from matplotlib import ticker as mticker
from matplotlib.ticker import MultipleLocator
rcParams.update({'font.size': 26}) # changed from 36 to 26
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
fig.subplots_adjust(top=0.985, bottom=0.15, left=0.19, right=0.97, hspace=0.05, wspace=0.2) # changed 0.17 to 0.19
#

# settings
folder = "/home/lukas/simulations/tdpot/project_thesis/performance/q8_10keV/"
file = os.path.join(folder, "performance_test.txt")
#

# get data
n = []
rk = []
vv = []

with open (file, "r") as f:
    lines = f.readlines()
    for line in lines:
        if not line.startswith("#"):
            n.append(float(line.split()[0]))
            rk.append(float(line.split()[1]))
            vv.append(float(line.split()[2]))
#

# plot

plt.plot(n, rk, color=t4iblue, label="Runge-Kutta method")
plt.plot(n, vv, color=t4imaroon, label ="Verlet integration")

plt.xlabel(r"Number of ions")
plt.ylabel(r"Duration of simulation (s)")

plt.legend()
plt.savefig(os.path.join(folder, "performance_test.jpg"), dpi=600)

    # from plot template
ax.xaxis.set_major_locator(MultipleLocator(10)) # changed from 5 to 10
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
#
