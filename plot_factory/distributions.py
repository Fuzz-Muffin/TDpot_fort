# plots the exit charge state and the angular deflection as histogram

import os

import matplotlib.pyplot as plt
import numpy as np

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
fig.subplots_adjust(top=0.985, bottom=0.16, left=0.19, right=0.97, hspace=0.05, wspace=0.2) # changed 0.15 to 0.16, changed 0.17 to 0.19
#

# settings
folder = "/home/lukas/simulations/tdpot/project_thesis/angles/prod"
target = "TLG"
energy = "100"
file_1 = f"{target}/Xe30_{energy}keV/{target}_Xe30_{energy}keV_g-0p0105-2p8-2-a-0/{target}_Xe30_{energy}keV_g-0p0105-2p8-2-a-0_out_2.txt"
file_2 = f"{target}/Xe30_{energy}keV/{target}_Xe30_{energy}keV_g-0p0105-2p8-2-a-10/{target}_Xe30_{energy}keV_g-0p0105-2p8-2-a-10_out_2.txt"
file_3 = f"{target}/Xe30_{energy}keV/{target}_Xe30_{energy}keV_g-0p0105-2p8-2-a-20/{target}_Xe30_{energy}keV_g-0p0105-2p8-2-a-20_out_2.txt"
#

# get data
qout = [[], [], []]
phi = [[], [], []]
psi = [[], [], []]
files = [file_1, file_2, file_3]
for counter, file in enumerate(files):
    with open(os.path.join(folder, file)) as f:
        for line in f:
            if not line.startswith(" #"):
                info = line.split()
                qout[counter].append(int(info[5]))
                phi[counter].append(np.arctan(float(info[7]))*180/np.pi)
                psi[counter].append(np.arctan(float(info[8]))*180/np.pi)
#

# plot
plt.hist(qout, histtype="step", bins=np.arange(0, 20, 1), color=[t4iblue, t4imaroon, t4ipetrol] , label=["0 deg", "10 deg", "20 deg"], stacked=False, fill=False, linewidth=2)
plt.xlabel(r"Charge state of the ion")
plt.ylabel(r"Number of ions")
plt.legend()
plt.savefig(os.path.join(folder, "hist_qout_" + target + "_" + energy + ".pdf"), dpi=600)

plt.clf()

    # setup
fig = plt.figure(figsize=(10.5, 7))
ax = fig.add_subplot(111)
plt.setp(ax.spines.values(), lw=3)
fig.subplots_adjust(top=0.985, bottom=0.16, left=0.19, right=0.97, hspace=0.05, wspace=0.2) # changed 0.15 to 0.16, changed 0.17 to 0.19
#

plt.hist(phi, histtype="step", bins=np.arange(0, 90, 2), color=[t4iblue, t4imaroon, t4ipetrol] , label=["0 deg", "10 deg", "20 deg"], stacked=False, fill=False, linewidth=2)
plt.xlabel(r"Angular deflection $\phi$ of the ion (degree)")
plt.ylabel(r"Number of ions")
plt.legend()
plt.savefig(os.path.join(folder, "hist_phi_" + target + "_" + energy + ".pdf"), dpi=600)

plt.clf()

    # setup
fig = plt.figure(figsize=(10.5, 7))
ax = fig.add_subplot(111)
plt.setp(ax.spines.values(), lw=3)
fig.subplots_adjust(top=0.985, bottom=0.16, left=0.19, right=0.97, hspace=0.05, wspace=0.2) # changed 0.15 to 0.16, changed 0.17 to 0.19
#

plt.hist(psi, histtype="step", bins=np.arange(0, 45, 2), color=[t4iblue, t4imaroon, t4ipetrol] , label=["0 deg", "10 deg", "20 deg"], stacked=False, fill=False, linewidth=2)
plt.xlabel(r"Angular deflection $\psi$ of the ion (degree)")
plt.ylabel(r"Number of ions")
plt.legend()
plt.savefig(os.path.join(folder, "hist_psi_" + target + "_" + energy + ".pdf"), dpi=600)

    # from plot template
ax.xaxis.set_major_locator(MultipleLocator(10)) # changed from 5 to 10
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
#
