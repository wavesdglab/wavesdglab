# WavesDGlab

Collection of MATLAB scripts developed to test and study high-order **continuous Galerkin (CG)** and **discontinuous Galerkin (DG)** finite element methods.

[![MATLAB](https://img.shields.io/badge/MATLAB-orange.svg?logo=mathworks&logoColor=white)](https://www.mathworks.com/products/matlab.html)
[![GMSH](https://img.shields.io/badge/mesh-GMSH-blue.svg)](https://gmsh.info/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Requirements

- MATLAB
- Gmsh

## Getting started

Before running the scripts, configure the path to the Gmsh installation:

1. Copy `setup_default.m` to `setup.m`.
2. Set the `directoryGmsh` variable in `setup.m` to the path of your Gmsh installation.
3. Run `setup.m` from the root directory of the repository.
4. Run one of the scripts located in the `scripts` directory.

## Repository structure

```
wavesdglab/
├── benchmarks/   # benchmark descriptions and data (meshes, reference solutions)
├── scripts/      # scripts to run experiments and simulations
├── src/          # source code
├── output/       # empty directory for output files (created at runtime)
├── tools/        # external resources and third-party contributions
└── README.md
```

## License

WavesDGlab is copyright © 2023–2026 CNRS, Inria, and ENSTA.

It is distributed under the terms of the **GNU General Public License, Version 3 (GPLv3)**.

**Exception:** The files in the `tools` directory are external contributions. See [`LICENSE.txt`](LICENSE.txt) for more information.

## Authors

- **Maintainer:** [Axel Modave](https://github.com/axmodave)
- **Contributors:** [Théophile Chaumont-Frelet](https://github.com/tchaumont), [Pierre Marchand](https://github.com/pierremarchand20), [Simone Pescuma](https://github.com/simonepescuma), [Timothée Raynaud](https://github.com/timotheeraynaud)

## References

**Hybridizable discontinuous Galerkin (HDG) methods**

- A. Modave, T. Chaumont-Frelet, *A hybridizable discontinuous Galerkin method with characteristic variables for Helmholtz problems*, **Journal of Computational Physics**, 493, 112459, 2023. [article](https://doi.org/10.1016/j.jcp.2023.112459) [preprint](https://hal.science/hal-03909368)
- S. Pescuma, G. Gabard, T. Chaumont-Frelet, A. Modave, *A hybridizable discontinuous Galerkin method with transmission variables for time-harmonic acoustic problems in heterogeneous media*, **Journal of Computational Physics**, 534, 114009, 2025. [article](https://doi.org/10.1016/j.jcp.2025.114009) [preprint](https://hal.science/hal-04821539)
- S. Pescuma, G. Gabard, T. Chaumont-Frelet, A. Modave, *A HDG method with transmission variables for time-harmonic wave propagation problems with constant coefficients*, 2026. [preprint](https://hal.science/hal-05654780)

**Continuous Galerkin (CG) methods**

- V. Dolean, P. Marchand, A. Modave, T. Raynaud, *Convergence analysis of GMRES applied to Helmholtz problems near resonances*, 2025. [preprint](https://hal.science/hal-05078654)