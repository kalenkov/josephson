# Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
# Licensed under the MIT License
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# path to .so (build/)
root = Path(__file__).resolve().parents[1]
sys.path.append(str(Path(__file__).resolve().parent))
sys.path.append(str(root / "build"))
sys.path.append(str(root / "build-cmake"))
import josephson

def print_config(stream, prefix, config, rt_info):
    print(f"{prefix}Kernel configuration:", file=stream)
    for key, value in config.items():
        print(f"{prefix}{key.replace('_','-')}={value}", file=stream)
    print(f"{prefix}Calculated parameters:", file=stream)
    print(f"{prefix}maximal grid energy={rt_info.grid_max_energy}", file=stream)
    print(f"{prefix}calculated step size={rt_info.step_size}", file=stream)
    print(f"{prefix}grid size={rt_info.grid_size}", file=stream)
    print("\n")

cfg = {
    "R_up": 0.1,
    "R_down": 0.1,
    "spin_mixing_angle": 0.0,
    "energy_cutoff": 20.0,
    "step_size": 0.01,
    "eV": 0.1,
    "T": 0.1,
    "energy_broadening": 1e-7,
    "order_parameter_cutoff_smoothing": 2.0,
    "show_progress": True,
}

jc = josephson.JosephsonCurrent(cfg)

print_config(sys.stdout, "# ", cfg, jc.runtime_info())

phi = np.linspace(0.0, 2.0*np.pi, 30)
J = jc.evaluate(phi)

print("# phi J")
for i in range(len(J)):
    print(phi[i], J[i])


f = interp1d(phi, J, kind='cubic')

phi_smooth = np.linspace(phi.min(), phi.max(), 500)
J_smooth = f(phi_smooth)

plt.plot(phi, J, 'o', markersize=3) 
plt.plot(phi_smooth, J_smooth, 'k-') 

plt.xlim(0, 2*np.pi)
xticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
xtick_labels = ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']
plt.xticks(xticks, xtick_labels)
plt.xlabel(r'$\phi$')
plt.ylabel('J')
plt.show()