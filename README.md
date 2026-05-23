# Josephson

## Description

Numerical framework for calculating the current-phase relation in nonequilibrium SIS symmetric Josephson junctions with a spin-active interface biased by a time-independent voltage.

## Features

- C++20 core library
- OpenMP parallelization
- Command-line interface
- Python bindings via pybind11

## Requirements

- C++20-compatible compiler
- CMake
- Boost
- Eigen3
- OpenMP
- Python 3
- pybind11

## Build instructions

Run the following commands to build the library and executables:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The project uses a standard CMake build system and produces all artifacts in the build directory. It has been tested on Fedora Linux with GCC and on Windows using the MSYS2 MinGW-w64 toolchain. The code is written in standard C++20 and is expected to be portable across modern Linux systems.

## Usage

The project can be used either via the command-line interface or through the Python API.

### Command-line interface

```text
josephson [OPTIONS]
```

| Option | Meaning | Default |
|---|---|---|
| `--R-up VALUE` | Spin-up reflection probability | `0.1` |
| `--R-down VALUE` | Spin-down reflection probability | `0.1` |
| `--spin-mixing-angle VALUE` | Spin-mixing angle | `0.0` |
| `--energy-cutoff VALUE` | Energy cutoff | `20` |
| `--step-size VALUE` | Energy grid step size | `0.001` |
| `--eV VALUE` | Applied voltage | `0.1` |
| `--T VALUE` | Temperature | `0.1` |
| `--energy-broadening VALUE` | Energy broadening parameter | `1e-7` |
| `--order-parameter-cutoff-smoothing VALUE` | Smoothing width of the order-parameter cutoff at high energies | `2.0` |
| `--phi-min VALUE` | Minimum phase | `0.0` |
| `--phi-max VALUE` | Maximum phase | `2π` |
| `--n-phi VALUE` | Number of phase points | `100` |
| `--output FILE` | Output filename | `output.txt` |
| `--show-progress` | Show computation progress | disabled |

### Example

The executable can be run as

```bash
./build/josephson \
    --R-up=0.1 \
    --R-down=0.1 \
    --spin-mixing-angle=0.0 \
    --eV=0.10 \
    --step-size=0.001 \
    --energy-broadening=1e-7 \
    --energy-cutoff=20 \
    --order-parameter-cutoff-smoothing=2 \
    --n-phi=3 \
    --T=0.1 \
    --phi-min=0.3 \
    --phi-max=1.5 \
    --show-progress
```

The same command-line interface is also available through the Python wrapper script:

```bash
python python/josephson_cli.py \
    --R-up=0.1 \
    --R-down=0.1 \
    --spin-mixing-angle=0.0 \
    --eV=0.10 \
    --step-size=0.001 \
    --energy-broadening=1e-7 \
    --energy-cutoff=20 \
    --order-parameter-cutoff-smoothing=2 \
    --n-phi=3 \
    --T=0.1 \
    --phi-min=0.3 \
    --phi-max=1.5 \
    --show-progress
```

The script demonstrates how to pass parameters to the `josephson` Python module, which internally uses the C++ core library for better performance.

A much simpler example is provided in `python/example_basic.py`.

Run

```bash
python python/example_basic.py
```

without any command-line arguments. The script demonstrates direct usage of the Python API, calculation of the current-phase relation, and plotting of the result. The parameters can be edited directly in the script.

## Parameter description

All energy scales (voltage \(eV\), temperature \(T\), energy cutoff, broadening, etc.) are measured in units of the superconducting order parameter magnitude \(|\Delta|\).

The superconducting order parameter is normalized such that \(|\Delta| 1\). This sets the natural energy scale of the problem.

### `R-up`, `R-down`

Reflection probabilities for spin-up and spin-down channels at the interface. These parameters vary from `0` (fully transparent interface) to `1` (fully reflecting interface).

They control the normal-state transparency of the junction for quasiparticles with different spin orientations.

### `spin-mixing-angle`

May take arbitrary values. Describes the difference between scattering phases acquired by spin-up and spin-down quasiparticles at the interface.

### `energy-cutoff`

The superconducting order parameter is set to zero for energies above this cutoff value. The parameter is measured in units of the superconducting order parameter.

This cutoff is introduced for numerical purposes. The calculated Josephson current becomes nearly independent of this parameter when it is chosen sufficiently large.

### `step-size`

The numerical algorithm discretizes the energy axis. This parameter controls the discretization resolution.

For reliable results, it should typically be smaller than both the temperature and the applied voltage.

The parameter is measured in units of the superconducting order parameter.

### `eV`

The parameter `eV` entering the theory, measured in units of the superconducting order parameter.

### `T`

Temperature measured in units of the superconducting order parameter.

### `energy-broadening`

Phenomenological quasiparticle damping parameter measured in units of the superconducting order parameter.

### `order-parameter-cutoff-smoothing`

An abrupt suppression of the order parameter at the cutoff energy may produce numerical instabilities. Therefore, the order parameter is smoothly reduced to zero within an energy interval specified by this parameter.

The parameter is measured in units of the superconducting order parameter.

### `phi-min`, `phi-max`, `n-phi`

The program evaluates the Josephson current at `n-phi` equally spaced phase values ranging from `phi-min` to `phi-max`.

### `output`

In addition to printing results to the console, save them to a file.

### `show-progress`

Some calculations may take a significant amount of time. This option enables printing of the estimated computation progress and remaining time.

## Additional notes

A brief summary of the theoretical model, normalization conventions, and equations implemented in the code can be found in

[Nonequilibrium Josephson Effect in SIS Junctions (PDF)](docs/notes.pdf)

## Citation

If this code is useful in your research, please cite

M.S. Kalenkov and A.D. Zaikin,  *Absence of fractional ac Josephson effect in superconducting junctions*, Phys. Rev. B **110**, 134521 (2024).

## Tests

To verify correctness of the build, run all tests with

```bash
ctest --test-dir build --output-on-failure
```

### C++ unit tests

Internal tests and correctness checks against reference points:

```bash
./build/josephson_test
```

### Python tests

Tests of the Python interface:

```bash
pytest -s python/test
```

### Notes

Tests may take a relatively long time depending on CPU performance.

## License

This project is distributed under the MIT License.
See the LICENSE file for details.

## Contact

For feedback, bug reports, or suggestions, please contact:

**Mikhail Kalenkov**  
Email: mikhail.kalenkov@gmail.com