# Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
# Licensed under the MIT License
import sys
import shlex
from pathlib import Path
import argparse
import numpy as np

# path to .so (build/)
root = Path(__file__).resolve().parents[1]
sys.path.append(str(Path(__file__).resolve().parent))
sys.path.append(str(root / "build"))
sys.path.append(str(root / "build-cmake"))
import josephson

def print_config(stream, prefix, config, rt_info, params):
    print(f"{prefix}Command line={params['cmd']}", file=stream)
    print(f"{prefix}Kernel configuration:", file=stream)
    for key, value in config.items():
        print(f"{prefix}{key.replace('_','-')}={value}", file=stream)
    print(f"{prefix}Calculated parameters:", file=stream)
    print(f"{prefix}maximal grid energy={rt_info.grid_max_energy}", file=stream)
    print(f"{prefix}calculated step size={rt_info.step_size}", file=stream)
    print(f"{prefix}grid size={rt_info.grid_size}", file=stream)
    print(f"{prefix}Program configuration:", file=stream)
    print(f"{prefix}phi-min={params['phi-min']}", file=stream)
    print(f"{prefix}phi-max={params['phi-max']}", file=stream)
    print(f"{prefix}n-phi={params['n-phi']}", file=stream)
    print(f"{prefix}filename={params['filename']}", file=stream)
    print("\n")

full_cmd = " ".join(shlex.quote(arg) for arg in sys.argv)

parser = argparse.ArgumentParser(description="Allowed options")
parser.add_argument("--R-up", type=float, default=0.1, help="Spin-up reflection probability")
parser.add_argument("--R-down", type=float, default=0.1, help="Spin-down reflection probability")
parser.add_argument("--spin-mixing-angle", type=float, default=0.0, help="Spin-mixing angle")
parser.add_argument("--energy-cutoff", type=float, default=20, help="Energy cutoff")
parser.add_argument("--step-size", type=float, default=0.001, help="Energy mesh step size")
parser.add_argument("--eV", type=float, default=0.1, help="Applied voltage")
parser.add_argument("--T", type=float, default=0.1, help="Temperature")
parser.add_argument("--energy-broadening", type=float, default=1e-07, help="Decay parameter")
parser.add_argument("--order-parameter-cutoff-smoothing", type=float, default=2, help="Smoothing width of the order-parameter cutoff at high energies")
parser.add_argument("--phi-min", type=float, default=0.0, help="Minimal phase")
parser.add_argument("--phi-max", type=float, default=6.283185307179586, help="Maximal phase")
parser.add_argument("--n-phi", type=int, default=100, help="Number of phase points")
parser.add_argument("-o", "--output", type=str, default="output.txt", help="Output filename")
parser.add_argument("--show-progress", action="store_true", help="Show computation progress")
args = parser.parse_args()
cfg = {
    "R_up": args.R_up,
    "R_down": args.R_down,
    "spin_mixing_angle": args.spin_mixing_angle,
    "energy_cutoff": args.energy_cutoff,
    "step_size": args.step_size,
    "eV": args.eV,
    "T": args.T,
    "energy_broadening": args.energy_broadening,
    "order_parameter_cutoff_smoothing": args.order_parameter_cutoff_smoothing,
    "show_progress": args.show_progress,
}

jc = josephson.JosephsonCurrent(cfg)

with open(args.output, "w") as f:
    streams = [sys.stdout, f]
    
    for s in streams:
        print_config(
            s,
            "# ",
            cfg,
            jc.runtime_info(),
            {
                "cmd": full_cmd,
                "phi-min": args.phi_min,
                "phi-max": args.phi_max,
                "n-phi": args.n_phi,
                "filename": args.output,
            },
        )

    phi = np.linspace(args.phi_min, args.phi_max, args.n_phi)
    J = jc.evaluate(phi)

    for s in streams:
        print("# phi J", file=s)
        for i in range(len(J)):
            print(phi[i], J[i], file=s)
