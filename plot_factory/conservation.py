import os

import matplotlib.pyplot as plt
import numpy as np

# from plot template
from matplotlib import rcParams
from matplotlib import ticker as mticker
from matplotlib.ticker import MultipleLocator
rcParams.update({'font.size': 9}) # changed from 36 to 9
rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})
rcParams['font.sans-serif'] = "Helvetica"
rcParams['font.family'] = "sans-serif"
rcParams['lines.linewidth'] = 1 # changed from 2 to 1
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
fig = plt.figure(figsize=(21, 7)) # changed from (10.5, 7) to (21, 7)
ax = fig.add_subplot(111)
plt.setp(ax.spines.values(), lw=3)
fig.subplots_adjust(top=0.985, bottom=0.15, left=0.19, right=0.97, hspace=0.05, wspace=0.2) # changed 0.17 to 0.19
#

# settings
folder = "/home/lukas/simulations/tdpot/project_thesis/performance/q8_10keV/"
method = "1"
ion_nr = "1"
#

# get data
energie_data = np.loadtxt(os.path.join(folder, "log_energies_" + method + "_" + ion_nr + ".txt"))
energie_data_T = energie_data.T

t = energie_data_T[0]
e_pot = energie_data_T[1]
e_kin = energie_data_T[2]
energy = energie_data_T[3]

electrons_data = np.loadtxt(os.path.join(folder, "log_e_number_" + method + "_" + ion_nr + ".txt"))
electrons_data_T = electrons_data.T

n_sta = electrons_data_T[2]
n_cap = electrons_data_T[3]
#

# plot
fig, host = plt.subplots(layout="constrained")

host.set_xlabel(r"Distance to target in z (a.u.)")
host.set_ylabel(r"$E_{pot}$ + $E_{kin}$ (eV)")
host.yaxis.label.set_color(t4imaroon)
ax2 = host.twinx()
ax2.set_ylabel(r"$E_{pot}$ (eV)")
ax2.yaxis.label.set_color(t4iblue)
ax3 = host.twinx()
ax3.set_ylabel(r"$E_{kin}$ (eV)")
ax3.yaxis.label.set_color(t4iorange)
ax4 = host.twinx()
ax4.set_ylabel(r"Captured electrons")
ax4.yaxis.label.set_color(t4ivio)
ax5 = host.twinx()
ax5.set_ylabel(r"Stabilized electrons")
ax5.yaxis.label.set_color(t4ipurp)

ax3.spines["right"].set_position(("outward", 40))
ax4.spines["right"].set_position(("outward", 100))
ax5.spines["right"].set_position(("outward", 160))

p1 = host.plot(t, energy, color=t4imaroon, linestyle=":", linewidth=2)
p2 = ax2.plot(t, e_pot, color=t4iblue)
p3 = ax3.plot(t, e_kin, color=t4iorange)
p4 = ax4.plot(t, n_cap, color=t4ivio, linestyle="--")
p5 = ax5.plot(t, n_sta, color=t4ipurp, linestyle="--")

plt.xlim(-15, 15)
plt.savefig(os.path.join(folder, "log_energies_" + method + "_" + ion_nr + ".jpg"), dpi=600)

    # from plot template
ax.xaxis.set_major_locator(MultipleLocator(10)) # changed from 5 to 10
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
#
