# Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
# Licensed under the MIT License

import sys
from pathlib import Path
import pytest
import csv

root = Path(__file__).resolve().parents[2]
sys.path.append(str(root / "build"))
datafile = root / "test" / "data.csv"

import josephson


def load_cases(fast_mode=False):
    cases = []
    with open(datafile, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if fast_mode and row["label"] != "fast":
                continue
            cases.append(
                {
                    "R_up": float(row["R_up"]),
                    "R_down": float(row["R_down"]),
                    "spin_mixing_angle": float(row["spin_mixing_angle"]),
                    "energy_cutoff": float(row["energy_cutoff"]),
                    "step_size": float(row["step_size"]),
                    "eV": float(row["eV"]),
                    "T": float(row["T"]),
                    "energy_broadening": float(row["energy_broadening"]),
                    "order_parameter_cutoff_smoothing": float(row["order_parameter_cutoff_smoothing"]),
                    "phi": float(row["phi"]),
                    "current": float(row["current"]),
                }
            )
    return cases


def pytest_generate_tests(metafunc):
    if "case" not in metafunc.fixturenames:
        return
    fast_mode = metafunc.config.getoption("--fast")
    metafunc.parametrize("case", load_cases(fast_mode=fast_mode))


def test_regression(case):
    cfg = {
        "R_up": case["R_up"],
        "R_down": case["R_down"],
        "spin_mixing_angle": case["spin_mixing_angle"],
        "energy_cutoff": case["energy_cutoff"],
        "step_size": case["step_size"],
        "eV": case["eV"],
        "T": case["T"],
        "energy_broadening": case["energy_broadening"],
        "order_parameter_cutoff_smoothing": case["order_parameter_cutoff_smoothing"],
        "show_progress": False,
    }
    jc = josephson.JosephsonCurrent(cfg)
    result = jc.evaluate(case["phi"])
    assert result == pytest.approx(
        case["current"],
        rel=1e-10,
    )
