# MD
This repository contains a **from-scratch molecular dynamics (MD) engine written in C++** to demonstrate MD algorithms and C++ coding proficiency. It numerically integrates [Newton’s equations](https://en.wikipedia.org/wiki/Newton%27s_laws_of_motion) for particles with bonded and non-bonded interactions using neighbor lists and periodic boundary conditions. Bonded terms include harmonic bonds, angles, and dihedrals; non-bonded terms include the [Lennard-Jones 12-6 potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) (for van der Waals forces) and [Coulombic electrostatics](https://en.wikipedia.org/wiki/Coulomb%27s_law). Long-range Coulomb forces are handled by the [Particle Mesh Ewald (PME)](https://en.wikipedia.org/wiki/Ewald_summation#Particle_mesh_Ewald_(PME)_method) method, which splits interactions into real-space and Fourier-space parts for rapid convergence.

It is a personal project built to demonstrate understanding of **molecular dynamics algorithms**, **numerical methods**, and **modern C++ design**, rather than a production-ready MD package.

The code implements core components found in real MD engines (bonded and non-bonded forces, neighbor lists, PME electrostatics, FFTs, periodic boundary conditions) with an emphasis on correctness, clarity, and performance-aware design.

> This project is **not intended for real scientific use**. It is a learning and portfolio project.

---

## Features

### Force Field Components

- **Bonded interactions**
  - Harmonic bond stretching
  - Angle bending
  - Dihedral (torsional) potentials

- **Non-bonded interactions**
  - Lennard-Jones 12-6 potential
  - Coulomb electrostatics

### Long-Range Electrostatics

- **Particle Mesh Ewald (PME)**
  - Real-space short-range electrostatics
  - Reciprocal-space electrostatics using FFTs
  - Charge spreading onto a mesh and force interpolation

- **FFT-based implementation**
  - Uses FFTW3 for fast Fourier transforms
  - O(N log N) scaling for long-range electrostatics

### Performance-Oriented Design

- **Neighbor / pair lists**
  - Verlet-style pair list construction
  - Reduces cost of short-range non-bonded interactions

- **Periodic boundary conditions**
  - Minimum image convention
  - Supports fully periodic simulation boxes

---
```
MD/
├── src/ # Core MD engine implementation
│ ├── AmberTopology/ # Amber Topology Reading
│ ├── Bonded/ # Bonded Interactions
│ ├── ChemicalEntity/ # Atoms, Molecules, Element
│ ├── NonBonded/ #NonBonded Interactions
│ ├── System/ #PBC and System Setup, Force Field Parameters and Box initialization
│ ├── Util/ Utilities for reading text, etc.
│
├── apps/ # Executables / drivers
│ ├── main.cpp #Main driver
│
├── CMakeLists.txt # CMake build configuration
└── README.md
```

## Installation

Prerequisites: Install CMake and the [FFTW3](https://www.fftw.org/) development library. It can be installed on Ubuntu-like systems using `sudo apt install libfftw3-dev`.

Clone and build:

```
git clone https://github.com/prateekbansal97/MD.git
cd MD
mkdir build && cd build
cmake ..
make
```

This produces an executable (e.g. main or MD) in the build directory.

Usage

Run the compiled program with a topology and coordinate file. For example, from the `build` directory:
```
./main ../CB1_apo_assym.prmtop ../CB1_apo_assym_fah_sub_100.rst7
```

This should produce an output of the form:
```
box: 66 66 66 90 90 90 
Box center: 20.8168 37.7142 26.1012 
Energies:
Bond:3617.4 
Angle: 2478.81 
CosineDihedral: 744.073 
Urey-Bradley: 34.9848 
Impropers: 183.599 
CMAP energy: -12.6229 
VDW: 12888.2 
EE: -123265 

Calculated Forces! Force on first three atoms (x, y, z directions): 
-1.56291 -3.23577 -5.86805 
-0.125699 -0.481944 -4.23094 
-2.59804 4.19055 -0.13331
```



(The included example files demonstrate one way to initialize a simulation.) The program will compute forces and energies based on the input structure. Note: This is a demonstration project; input/output formats and user interface are minimal. Users should examine the source in apps/ or src/ to understand how the simulation is driven.
Author & License

This code was solely written by Prateek Bansal as a personal project to showcase MD implementation skills. It is released under the MIT License (see the LICENSE file for details).

