# Waves in Granular Materials

Numerical simulations of wave propagation in two-dimensional granular materials using a strongly nonlinear continuous model.

This repository contains MATLAB implementations of numerical methods for solving the governing equations of granular media, based on the three-contact constitutive model developed in Cominelli et al.

The project focuses on the propagation of nonlinear waves in granular materials, accounting for particle interactions, deformation, stress, gravity, and dissipation.

## Scientific Background

The code implements a continuous model describing the macroscopic behavior of two-dimensional granular materials.

The model considers a granular assembly in which the mechanical energy originates from three contact directions. The governing equations describe the evolution of:

- **Density** (\(\rho\))
- **Velocity field** (\(u,v\))
- **Deformation gradient** (\(F\))
- **Stress tensor** (\(\sigma\))

The constitutive behavior is based on Hertzian contact mechanics.

The numerical formulation uses operator splitting to separate the governing equations into two subsystems:

1. A subsystem describing mass and momentum conservation, solved using Bale's numerical scheme.
2. A subsystem describing the evolution of the deformation-related variables, solved using a Rusanov-type scheme or related numerical methods.

The implementation also includes adaptive time stepping controlled by a CFL condition.

## Repository Structure

The main files are organized as follows:

| File | Description |
|---|---|
| `main_Comi_3_contacts_splitting.m` | Main simulation script. Defines the physical parameters, initial conditions, boundary conditions, numerical schemes, and executes the simulation. |
| `solver.m` | General numerical solver and time-integration loop. |
| `set_scheme.m` | Configuration of the numerical schemes. |
| `Bale.m` | Implementation of Bale's numerical scheme. |
| `Rusanov.m` | Rusanov numerical flux implementation. |
| `MUSCL_H.m` | MUSCL-Hancock reconstruction and numerical scheme. |
| `Runge_Kutta.m` | Runge-Kutta time integration routines. |
| `app_BCs.m` | Application of boundary conditions. |
| `A_fun.m`, `B_fun.m` | Matrices associated with the governing system. |
| `fun_cont_hertz_homo.m` | Hertzian constitutive model for homogeneous granular materials. |
| `c_fun_Comi_3_contacts.m` | Characteristic wave speeds for the three-contact model. |
| `eigs_x_3_contacts.m` | Eigenvalue and eigenvector calculations. |
| `estrai.m` | Extraction of model variables and auxiliary quantities. |
| `cosa_plot_3_contacts.m` | Visualization of simulation results. |
| `fig_bella.m` | Figure generation and visualization utilities. |
| `save_GIF.m` | GIF generation from simulation frames. |

The repository also contains `.mat` files with simulation data, including configurations associated with the Melin granular packing.

## Requirements

- MATLAB
- A MATLAB installation capable of running the scripts and numerical routines included in this repository.

No external package manager or build system is required.

The code has been implemented using MATLAB R2025a.

## Installation

Clone the repository:

```bash
git clone https://github.com/SebaComi/Waves_granular_material.git
```

Open the project directory in MATLAB or add it to the MATLAB path.

## Running a Simulation

The main entry point is:

```matlab
main_Comi_3_contacts_splitting
```

The main script initializes the physical model, constructs the numerical fluxes, defines the initial and boundary conditions, and calls the general solver.

### Simulation Configuration

The main simulation parameters can be modified in `main_Comi_3_contacts_splitting.m`.

Some important parameters include:

| Parameter | Description |
|---|---|
| `DATA.T` | Final simulation time. |
| `DATA.CFL` | CFL number controlling the adaptive time step. |
| `DATA.Nx`, `DATA.Ny` | Number of computational grid points, including ghost cells. |
| `DATA.R` | Particle radius. |
| `DATA.E_p` | Young's modulus of the particle material. |
| `DATA.nu_p` | Poisson's ratio. |
| `DATA.precompr` | Initial precompression parameter. |
| `DATA.fric` | Dissipation function. |
| `DATA.vv` | Imposed boundary velocity. |
| `DATA.scheme_type` | Numerical schemes used for the split subsystems. |

The default configuration uses:

```matlab
DATA.CFL = 0.6;
DATA.Nx = 50 + 4;
DATA.precompr = 1 + 2.15e-4;
DATA.T = 0.05;

DATA.scheme_type = {'Bale','MUSCL_H2D'};
```

These values can be adjusted to investigate different physical configurations and numerical behaviors.

### Initial Conditions and Data Import

By default, the main script initializes a homogeneous granular configuration.

To initialize a simulation from previously saved data, modify the `importa` flag:

```matlab
importa = true;
```

When enabled, the script loads a `.mat` file containing previously computed simulation data. Make sure the selected file exists in the working directory or update its path accordingly.

### Visualization

The simulation includes routines for visualizing physical quantities such as:

- Density.
- Velocity components.
- Stress tensor components.
- Deformation gradient components.
- Kinetic and potential energy.
- Total mass and energy evolution.

The `plots` and `salva` flags control visualization and GIF generation in the main script.

## Numerical Methods

The numerical framework is based on a splitting of the governing equations.

The conservative state vector is:

\[
W =
\begin{bmatrix}
\rho &
\rho u &
\rho v &
\rho F_{11} &
\rho F_{12} &
\rho F_{21} &
\rho F_{22}
\end{bmatrix}^{T}
\]

The first subsystem evolves the density and momentum variables, while the second subsystem evolves the deformation-related variables.

The implementation combines:

- Bale's scheme for the first subsystem.
- Rusanov-type numerical fluxes and MUSCL-Hancock methods for the second subsystem.
- Runge-Kutta integration routines.
- CFL-controlled adaptive time stepping.
- Wall boundary conditions.

The solver is designed to investigate nonlinear wave propagation and the behavior of the continuous granular model.

## Scientific Reference

This repository is associated with the following scientific work:

**Cominelli et al. — Strongly non-linear continuous model for wave propagation in 2D granular media: derivation and hyperbolicity analysis.**

The numerical implementation should be considered together with the theoretical formulation presented in this paper.

If you use this code in academic work, please cite the associated publication.

## License

This project is distributed under the **GNU General Public License v3.0 (GPL-3.0)**.

See the [`LICENSE`](LICENSE) file for the complete license terms.

## Author

**SebaComi**

GitHub: [@SebaComi](https://github.com/SebaComi)

---

*This repository is intended for scientific research and numerical investigation of nonlinear wave propagation in granular materials.*
