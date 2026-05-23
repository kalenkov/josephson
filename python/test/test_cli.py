# Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
# Licensed under the MIT License
from pathlib import Path
import subprocess
import sys
import tempfile
import math

import pytest


ROOT = Path(__file__).resolve().parents[2]
CLI_DIR = ROOT / "python" / "test" / "cli"
CPP_EXECUTABLE = ROOT / "build" / "josephson"
PYTHON_SCRIPT = ROOT / "python" / "josephson_cli.py"


FLOAT_KEYS = {
    "R-up",
    "R-down",
    "spin-mixing-angle",
    "energy-cutoff",
    "step-size",
    "eV",
    "T",
    "energy-broadening",
    "order-parameter-cutoff-smoothing",
    "maximal grid energy",
    "calculated step size",
}


BOOL_KEYS = {
    "show-progress",
}

OUTPUT_FILE_KEY = {"filename"}


CASES = [
    "case1",
]


def normalize_bool(value):
    value = value.strip().lower()
    if value == "true":
        return "true"
    if value == "false":
        return "false"
    raise ValueError(f"Invalid boolean value: {value}")


def normalize_line(line):
    line = line.strip()
    if "=" not in line:
        return line
    left, right = line.split("=", 1)
    left = left.strip()
    right = right.strip()
    key = left.lstrip("#").strip()
    if key in FLOAT_KEYS:
        return [left, float(right)]
    if key in OUTPUT_FILE_KEY:
        right = "output.txt"
    if key in BOOL_KEYS:
        right = normalize_bool(right)
    return f"{left} = {right}"


def parse_output(filename):

    header = []

    phi = []
    current = []

    in_data = False

    with open(filename) as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line == "# phi J":
                in_data = True
                continue
            if in_data:
                x, y = map(float, line.split())
                phi.append(x)
                current.append(y)
                continue
            if line.startswith("# Command line"):
                continue
            header.append(normalize_line(line))
    return header, phi, current


def compare_outputs(ref_file, test_file, atol=1e-12):
    ref_header, ref_phi, ref_current = parse_output(ref_file)
    test_header, test_phi, test_current = parse_output(test_file)
    assert len(ref_header) == len(test_header)
    for i, (ref, test) in enumerate(zip(ref_header, test_header)):
        if isinstance(ref, list):
            assert isinstance(test, list)
            assert ref[0] == test[0], \
                f"Header line {i}: {ref[0]} != {test[0]}"
            assert math.isclose(
                ref[1],
                test[1],
                abs_tol=atol,
            ), \
                f"Header float line {i}: {ref[1]} != {test[1]}"
        else:
            assert ref == test, \
                f"Header line {i}: {ref} != {test}"
    assert len(ref_phi) == len(test_phi)
    for i, (a, b) in enumerate(zip(ref_phi, test_phi)):
        assert math.isclose(a, b, abs_tol=atol), \
            f"phi[{i}]: {a} != {b}"
    assert len(ref_current) == len(test_current)
    for i, (a, b) in enumerate(zip(ref_current, test_current)):
        assert math.isclose(a, b, abs_tol=atol), \
            f"J[{i}]: {a} != {b}"


@pytest.mark.parametrize("case", CASES)
def test_cpp_cli(case):
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = Path(tmpdir) / "output.txt"
        args = (
            (CLI_DIR / f"{case}.args")
            .read_text()
            .split()
        )
        args += ["--output", str(outfile)]
        subprocess.run(
            [str(CPP_EXECUTABLE)] + args,
            check=True,
        )
        compare_outputs(
            CLI_DIR / f"{case}.ref",
            outfile,
        )


@pytest.mark.parametrize("case", CASES)
def test_python_cli(case):
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = Path(tmpdir) / "output.txt"
        args = (
            (CLI_DIR / f"{case}.args")
            .read_text()
            .split()
        )
        args += ["--output", str(outfile)]
        subprocess.run(
            [
                sys.executable,
                str(PYTHON_SCRIPT),
            ] + args,
            check=True,
        )
        compare_outputs(
            CLI_DIR / f"{case}.ref",
            outfile,
        )