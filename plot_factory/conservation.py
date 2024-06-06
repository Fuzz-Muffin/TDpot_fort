import os

import matplotlib.pyplot as plt
import numpy as np

folder = "/home/lukas/simulations/tdpot/slg/xe/q8/ekin167/"
method = "rk"

energie_data = np.loadtxt(os.path.join(folder, "log_energies_" + method + ".txt"))
energie_data_T = energie_data.T

t = energie_data_T[0]
e_pot = energie_data_T[1]
e_kin = energie_data_T[2]
energy = energie_data_T[3]

electrons_data = np.loadtxt(os.path.join(folder, "log_e_number_" + method + ".txt"))
electrons_data_T = electrons_data.T

n_sta = electrons_data_T[2]
n_cap = electrons_data_T[3]

if False:
    e_pot /= np.max(np.abs(e_pot))
    e_kin /= np.max(np.abs(e_kin))
    energy /= np.max(np.abs(energy))

if False:
    e_pot /= e_pot[0]
    e_kin /= e_kin[0]
    energy /= energy[0]

fig, host = plt.subplots(layout="constrained")

host.set_xlabel("count")
host.set_ylabel("E_pot + E_kin (eV)")

ax2 = host.twinx()
ax2.set_ylabel("E_pot (eV)")

ax3 = host.twinx()
ax3.set_ylabel("E_kin (eV)")

ax4 = host.twinx()
ax4.set_ylabel("e sta")

ax5 = host.twinx()
ax5.set_ylabel("e cap")

p1 = host.plot(t, energy, label="E_pot + E_kin", color="blue")
p2 = ax2.plot(t, e_pot, label="E_pot", color="orange", linestyle="--")
p3 = ax3.plot(t, e_kin, label="E_kin", color="green", linestyle="--")
ax3.spines["right"].set_position(("outward", 50))

p4 = ax4.plot(t, n_sta, label="e sta", color="red")
ax4.spines["right"].set_position(("outward", 120))
p5 = ax5.plot(t, n_cap, label="e cap", color="yellow")
ax5.spines["right"].set_position(("outward", 170))

plt.xlim(4000, 7001)
host.legend(handles=p1+p2+p3+p4+p5, loc="lower right")

plt.savefig(os.path.join(folder, "log_energies_" + method + ".jpg"), dpi=600)
