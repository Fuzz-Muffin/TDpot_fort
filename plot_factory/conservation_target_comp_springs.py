# compares the conservation of kinetic energy in the target in respect of springs

import os

import matplotlib.pyplot as plt
import numpy as np

# from plot template
from matplotlib import rcParams
from matplotlib import ticker as mticker
from matplotlib.ticker import MultipleLocator
rcParams.update({'font.size': 20}) # changed from 36 to 20
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
fig.subplots_adjust(top=0.985, bottom=0.15, left=0.1, right=0.97, hspace=0.05, wspace=0.2) # changed 0.17 to 0.1
#

# settings
folder_1 = "/home/lukas/simulations/tdpot/project_thesis/springs/q30_100keV/tlg/"
folder_2 = "/home/lukas/simulations/tdpot/project_thesis/springs/q30_100keV/tlg_no_springs/"
integrator = "2"
ion_nr = "2"
#

# get data
energie_data_1 = np.loadtxt(os.path.join(folder_1, "log_energies_" + integrator + "_" + ion_nr + ".txt"))
energie_data_1_T = energie_data_1.T

t_1 = energie_data_1_T[0]
e_kin_1 = energie_data_1_T[5]

energie_data_2 = np.loadtxt(os.path.join(folder_2, "log_energies_" + integrator + "_" + ion_nr + ".txt"))
energie_data_2_T = energie_data_2.T

t_2 = energie_data_2_T[0]
e_kin_2 = energie_data_2_T[5]
#

# plot
plt.plot(t_1, e_kin_1, color=t4iblue, label="Harmonic bonds")
plt.plot(t_2, e_kin_2, color=t4imaroon, label ="No harmonic bonds")

plt.xlabel(r"Distance to target in z (a.u.)")
plt.ylabel(r"$E_{kin}$ (eV)")

plt.legend()
plt.savefig(os.path.join(folder_2, "log_energies_bonds_" + integrator + "_" + ion_nr + ".jpg"), dpi=600)

    # from plot template
ax.xaxis.set_major_locator(MultipleLocator(10)) # changed from 5 to 10
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
#
