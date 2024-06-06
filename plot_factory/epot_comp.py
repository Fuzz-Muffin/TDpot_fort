import os

import matplotlib.pyplot as plt
import numpy as np

folder = "/home/lukas/simulations/tdpot/slg/xe/q8/ekin167/"

energie_data_rk = np.loadtxt(os.path.join(folder, "log_energies_rk.txt"))
energie_data_rk_T = energie_data_rk.T

energie_data_vv = np.loadtxt(os.path.join(folder, "log_energies_vv.txt"))
energie_data_vv_T = energie_data_vv.T

t = energie_data_rk_T[0]
energy_rk = energie_data_rk_T[1]
energy_vv = energie_data_vv_T[1]

if False:
    energy_rk /= np.max(np.abs(energy_rk))
    energy_vv /= np.max(np.abs(energy_vv))

if False:
    energy_rk /= energy_rk[0]
    energy_vv /= energy_vv[0]

plt.xlabel("count")
plt.ylabel("E_pot (eV)")

plt.plot(t, energy_rk, label="rk")
plt.plot(t, energy_vv, label="vv", linestyle="--")

plt.xlim(4000, 7001)
plt.legend()

plt.savefig(os.path.join(folder, "epot_comp.jpg"), dpi=300)
