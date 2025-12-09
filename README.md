# Reproduction material for the paper 'An unfitted divergence-free higher order finite element method for the Stokes problem' by M. Neilan, M. Olshanskii and H. von Wahl

[![DOI](https://zenodo.org/badge/1106042778.svg)](https://doi.org/10.5281/zenodo.17871850)

This repository contains the reproduction material for the paper 'An unfitted divergence-free higher order finite element method for the Stokes problem' by M. Neilan, M. Olshanskii and H. von Wahl. This includes a cpp add-on module to NGSolve [1] implementing the finite element space used in the paper. This add-on is based on [2-4]. The material in [4] is used under the MIT license, see `module/src/LICENSE-NGS-myfe`.

The repository contains the following files and directories:

| File/Directory | Description 
| --- | --- 
| `module/` | Contains the cpp add-on module implementing the FE space.
| `module/src/*` | Contains the finite element space implementation.
| `module/demo/test_spaces_stokes.py` | Contains a simple demo of the finite element space for a *fitted* Stokes problem.
| `run_examples.bash` | Bash script to run the examples presented in the paper.
| `stokes.py` | Main implementation of the unfitted method.
| `convergence.py` | Wrapper for convergence studies with given parameters.
| `condition.py` | Wrapper for condition number tests.
| `estcond.py` | Function to approximate the condition number of a linear system.
| `example0_data.py` | Problem data for example on circular geometry.
| `example1_data.py` | Problem data for Example 1.
| `example2_data.py` | Problem data for Example 2.
| `results/*` | Raw results presented in the paper.
| `README.md` | This file.
| `LICENSE` | License file.

### Installation
To run the python scripts locally, a compatible combination of `Netgen/NGSolve` and `ngsxfem` are required, together with the add-on module provided here. `Netgen/NGSolve` and `ngsxfem` can be installed by building from sources or precompiled pip wheels. For detailed installation instructions, we refer to the installation guidelines of [NGSolve](https://docu.ngsolve.org/latest/install/install_sources.html), and [ngsxfem](https://github.com/ngsxfem/ngsxfem/blob/release/INSTALLATION.md). Our numerical results are realised using the following versions:

| Package | git commit
|-|-|
| NGSolve | `a28640189a688b9e98ece17930fe3cd90688fcdd`
| ngsxfem | `86b395d2a5173faedee7b67edc4db569aef3358f`

To install the add-on module, we recommend that you build `NETGEN/NGSolve` and `ngsxfem` from their sources to ensure that all pre-requisits are available. The installation of the add-on module can then be compiled and installed using `cmake` as follows:
```
cd module
mkdir build
cd build
cmake ..
make -j4 install
```

### References
[1] https://ngsolve.org

[2] https://github.com/NGSolve/ngsolve-addon-template

[3] https://github.com/NGSolve/mylittlengsolve

[4] https://github.com/TUWien-ASC/NGS-myfe