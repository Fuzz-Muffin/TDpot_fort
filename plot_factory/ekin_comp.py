import os

import matplotlib.pyplot as plt
import numpy as np

folder = "/home/lukas/simulations/tdpot/slg/xe/q8/ekin167/"

energy_data = np.loadtxt(os.path.join(folder, "log_energies.txt"))
energy_data_T = energy_data.T

t = energy_data_T[0]
e_kin = energy_data_T[2]
energy_ion = energy_data_T[4]
energy_target = energy_data_T[5]

fig, host = plt.subplots(layout="constrained")

host.set_xlabel("count")
host.set_ylabel("E_kin (eV)")

ax2 = host.twinx()
ax2.set_ylabel("E_kin ion (eV)")

ax3 = host.twinx()
ax3.set_ylabel("E_kin target (eV)")

p1 = host.plot(t, e_kin, label="E_kin", color="blue")
p2 = ax2.plot(t, energy_ion, label="E_kin ion", color="orange", linestyle="--")
p3 = ax3.plot(t, energy_target, label="E_kin target", color="green", linestyle="--")
ax3.spines["right"].set_position(("outward", 60))

plt.xlim(4000, 7001)
host.legend(handles=p1+p2+p3, loc="center left")

plt.savefig(os.path.join(folder, "ekin_comp.jpg"), dpi=300)
